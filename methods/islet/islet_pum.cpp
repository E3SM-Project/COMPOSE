#include "islet_pum.hpp"

#include <omp.h>

#include "islet_util.hpp"
#include "islet_isl.hpp"

namespace pum {

static void init_xb (const Int ne, const Real& perturb, Real* const xb) {
  for (Int i = 0; i <= ne; ++i)
    xb[i] = 2*(Real(i)/ne) - 1;
  for (Int i = 1; i < ne; ++i)
    xb[i] += (perturb/ne)*(2*islet::urand() - 1);
}

static void fill_xgrid (const Int ne, const Real* const xb,
                        const Int np, const Real* const xnodes,
                        Real* const x) {
  for (Int i = 0; i < ne; ++i)
    for (Int j = 0; j < np; ++j) {
      const auto alpha = (1 + xnodes[j])/2;
      x[i*(np - 1) + j] = (1 - alpha)*xb[i] + alpha*xb[i+1];
    }
}

static void prep_for_eval (const Int ne, const Real* const xb, const Int np,
                           const Real& x, Int& ie, Real& xr, Real* const v) {
  for (ie = 0; ie < ne; ++ie)
    if (x >= xb[ie] && (ie == ne-1 || x < xb[ie+1]))
      break;
  assert(ie >= 0 && ie < ne);
  xr = 2*(x - xb[ie]) / (xb[ie+1] - xb[ie]) - 1;
  std::fill(v, v + ne*(np-1) + 1, 0);
}

Options::Options ()
  : threaded(false),
    ne(3), ntrial(111), mec_ne(11),
    perturb(0.05)
{}

PerturbedUniformMeshMetric
::PerturbedUniformMeshMetric (const InterpMethod::Ptr& im, const Options opts_)
  : opts(opts_), base_im(im), mec(opts.threaded)
{ init(); }

PerturbedUniformMeshMetric
::PerturbedUniformMeshMetric (const UserInterpMethod::Ptr& uim, const Options opts_)
  : opts(opts_), mec(opts.threaded)
{
  base_im = std::make_shared<InterpMethod>(uim);
  init();
}

void PerturbedUniformMeshMetric::reset_opts (const Options& o) {
  opts = o;
  init();
}

void PerturbedUniformMeshMetric::init () {
  xbs.resize(opts.threaded ? omp_get_max_threads() : 1);
  for (auto& xb : xbs) xb.resize(opts.ne+1);
  xnodess.resize(opts.threaded ? omp_get_max_threads() : 1);
  for (auto& xnodes : xnodess) xnodes.resize(get_np());
}

void PerturbedUniformMeshMetric::eval(const Real& x, Real* const v) {
  Int ie;
  Real xr;
  const Int tid = opts.threaded ? omp_get_thread_num() : 0;
  prep_for_eval(opts.ne, xbs[tid].data(), base_im->np,
                x, ie, xr, v);
  op_eval(*base_im, xr, v + ie*(base_im->np-1));
}

const Real* PerturbedUniformMeshMetric::get_xnodes() const {
  const Int tid = opts.threaded ? omp_get_thread_num() : 0;
  return xnodess[tid].data();
}

Int PerturbedUniformMeshMetric::get_np() const {
  return (xbs[0].size() - 1)*(base_im->np - 1) + 1;
}

Real PerturbedUniformMeshMetric::run (const Real stop_if_above,
                                      const bool one_elem_hop_only) {
  assert(base_im->np > 0);
  std::vector<Real> dxs; {
    dxs.push_back(1.0/opts.ne);
    if ( ! one_elem_hop_only) {
      for (const auto dx : {2.0/opts.ne, 0.5/opts.ne, 0.5})
        dxs.push_back(dx);
      const Real* xnodes = base_im->get_xnodes();
      for (int i = 0; i < base_im->np; ++i) {
        const Real xn = xnodes[i];
        if (xn == 1 || xn == -1) continue;
        dxs.push_back((0.5*(1 - xn))/opts.ne);
      }
    }
  }
  const auto run1 = [&] () {
    const int tid = opts.threaded ? omp_get_thread_num() : 0;
    auto& xb = xbs[tid];
    auto& xnodes = xnodess[tid];
    init_xb(opts.ne, opts.perturb, xb.data());
    fill_xgrid(opts.ne, xb.data(), base_im->np, base_im->get_xnodes(),
               xnodes.data());
    MaxEigComputer mec(false);
    Real mea, max_mea = 0;
    for (const auto dx : dxs) {
      mec.compute(this, dx, opts.mec_ne, &mea);
      max_mea = std::max(max_mea, mea);
      if (max_mea > 1 + stop_if_above) break;
    }
    return max_mea;
  };
  Real max_mea = 0;
  if (opts.threaded) {
#   pragma omp parallel for
    for (Int trial = 0; trial < opts.ntrial; ++trial) {
      if (max_mea > 1 + stop_if_above) continue;
      const auto mea = run1();
      if (mea > max_mea) {
#       pragma omp critical (PerturbedUniformMeshMetric_run)
        max_mea = std::max(max_mea, mea);
#       pragma omp flush
      }
    }
  } else {
    for (Int trial = 0; trial < opts.ntrial; ++trial) {
      const auto mea = run1();
      max_mea = std::max(max_mea, mea);
      if (max_mea > 1 + stop_if_above) break;
    }    
  }
  return max_mea - 1;
}

void PerturbedUniformMeshMetric
::sweep_and_collect_amplitudes (
  const Int npts, const Real threshold, std::map<Real,Real>& dx2meam1,
  const bool verbose)
{
  std::vector<Real> dxs; {
    for (const auto dx : {1.0/opts.ne, 0.5/opts.ne})
      dxs.push_back(dx);
    const Real* xnodes = base_im->get_xnodes();
    for (Int i = 0; i < base_im->np; ++i) {
      const Real xn = xnodes[i];
      if (xn == 1 || xn == -1) continue;
      dxs.push_back((0.5*(1 - xn))/opts.ne);
    }
    for (Int i = 1; i < npts; ++i)
      dxs.push_back(Real(i)/(opts.ne*npts));
  }
  const Int ndx = dxs.size();
# pragma omp parallel for schedule(static,1)
  for (Int trial = 0; trial < opts.ntrial; ++trial) {
    const int tid = omp_get_thread_num();
    auto& xb = xbs[tid];
    auto& xnodes = xnodess[tid];
    init_xb(opts.ne, opts.perturb, xb.data());
    fill_xgrid(opts.ne, xb.data(), base_im->np, base_im->get_xnodes(), xnodes.data());
    for (Int i = 0; i < ndx; ++i) {
      MaxEigComputer mec(false /* threaded */);
      Real mea;
      mec.compute(this, dxs[i], opts.mec_ne, &mea);
      if (mea >= 1 + threshold) {
#       pragma omp critical
        {
          const Real dx = dxs[i]*opts.ne;
          bool insert = false;
          if (dx2meam1.find(dx) != dx2meam1.end()) {
            const Real prev = dx2meam1[dx];
            if (mea > 1 + prev) {
              dx2meam1[dx] = mea-1;
              insert = true;
            }
          } else {
            dx2meam1[dx] = mea-1;
            insert = true;
          }
          if (insert && verbose) printf("dx %1.16e meam1 %1.16e\n", dx, mea-1);
        }
      }
    }
  }
}

void demo () {
  for (const Int np : {4, 5, 6, 7, 8, 9, 10, 11, 12, 13}) {
    Real pum_meam1_m0 = 0;
    for (Int method = 0; method < 2; ++method) {
      if (method == 1 && (np == 9 || np == 11)) continue;
      const auto oim = std::make_shared<islet::OperatorInterpMethod>(
        np, islet::Operator::create(method == 0 ?
                                    islet::Operator::gll_offset_nodal_subset :
                                    islet::Operator::xnodal));
      const auto im = std::make_shared<InterpMethod>(oim);
      Real pum_meam1; {
        pum::Options o;
        o.threaded = true;
        o.ne = 4;
        o.ntrial = 71;
        o.mec_ne = 3;
        o.perturb = 0.01;
        pum::PerturbedUniformMeshMetric pum(im, o);
        pum_meam1 = pum.run();
      }
      if (method == 0) pum_meam1_m0 = pum_meam1;
      Real meam1 = -1; {
        // Check meam1 to be sure we transcribed the new methods correctly.
        MaxEigComputer mec;
        const Int ns = 4111;
        meam1 = mec.run(np, ns, ns, 1e-14, true, im->uim);
      }
      printf("np %2d method %d meam1 %10.3e pum_meam1 %10.3e",
             np, method, meam1, pum_meam1);
      if (method == 0) printf("\n");
      else printf(" better %6.2f\n", pum_meam1_m0/pum_meam1);
    }
  }
}

} // namespace pum
