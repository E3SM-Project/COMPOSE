#include <set>

#include "slmm_mesh.hpp"
#include "slmm_accum.hpp"

#include "slmmir_remapper.hpp"
#include "slmmir_util.hpp"
#include "slmmir_p_refine.hpp"
#include "slmmir_d2c.hpp"

struct SearchFunctor {
  static constexpr Real pad = 1e-5;
  const Mesh& m;
  const Real* const p;
  Real bb[6];
  int ci_hit;

  SearchFunctor (const Mesh& m_, const Real p_[3])
    : m(m_), p(p_)
  {
    for (int j = 0; j < 3; ++j) bb[j  ] = p[j] - pad;
    for (int j = 0; j < 3; ++j) bb[j+3] = p[j] + pad;
    ci_hit = -1;
  }

  void operator() (const Int ci) {
    // Asssure there's always a hit in case numerics makes 'inside' return false
    // for all potential hits.
    if (ci_hit == -1) ci_hit = ci;
    const auto e = slice(m.geo_c2n, ci);
    const auto en = slice(m.geo_c2nml, ci);
    for (int i = 0; i < 4; ++i) {
      const Real* const a = slice(m.geo_p, e[i]);
      const Real* const nml = slice(m.geo_nml, en[i]);
      if ( ! siqk::SphereGeometry::inside(p, a, nml)) return;
    }
    ci_hit = ci;
  }
};

namespace sphere2ref {
// Better impl of sphere::siqk solver. Should move this code to that.
using namespace siqk;
using sqr::Info;

KOKKOS_INLINE_FUNCTION
void solve_Jxr (Real J[6], const Real r[3], Real dx[2]) {
  // QR factorization: J -> J [n1 a; 0 n2].
  const Real n1 = std::sqrt(SphereGeometry::norm2(J));
  SphereGeometry::scale(1/n1, J);
  const Real a = SphereGeometry::dot(J, J+3);
  SphereGeometry::axpy(-a, J, J+3);
  const Real n2 = std::sqrt(SphereGeometry::norm2(J+3));
  SphereGeometry::scale(1/n2, J+3);
  // r -> Q' r.
  Real Qtr[2] = {0};
  for (Int j = 0; j < 2; ++j) {
    const Real* const Jj = J + 3*j;
    for (Int i = 0; i < 3; ++i)
      Qtr[j] += Jj[i]*r[i];
  }
  // dx = R \ (Q' r).
  dx[1] = (Qtr[1] / n2);
  dx[0] = ((Qtr[0] - a*dx[1]) / n1);
}

KOKKOS_INLINE_FUNCTION
void calc_Jacobian (const Real c[12], const Real& a, const Real& b,
                    Real s_ab[6], Real s[3]) {
  static constexpr Real qtr = 0.25;
  {
    Real q[4];
    q[0] = qtr*(1-a)*(1-b); q[1] = qtr*(1+a)*(1-b);
    q[2] = qtr*(1+a)*(1+b); q[3] = qtr*(1-a)*(1+b);
    for (Int d = 0; d < 3; ++d) {
      const Real* const cd = c + 4*d;
      Real accum = 0;
      for (Int i = 0; i < 4; ++i) accum += cd[i]*q[i];
      s[d] = accum;
    }
  }
  const Real r2 = SphereGeometry::dot(s, s), r = std::sqrt(r2);
  for (Int j = 0; j < 2; ++j) {
    Real q_ab_j[4];
    if (j == 0) {
      q_ab_j[0] = -qtr*(1-b); q_ab_j[1] =  qtr*(1-b);
      q_ab_j[2] =  qtr*(1+b); q_ab_j[3] = -qtr*(1+b);
    } else {
      q_ab_j[0] = -qtr*(1-a); q_ab_j[1] = -qtr*(1+a);
      q_ab_j[2] =  qtr*(1+a); q_ab_j[3] =  qtr*(1-a);
    }
    Real* const s_ab_j = s_ab + 3*j;
    for (Int d = 0; d < 3; ++d) {
      const Real* const cd = c + 4*d;
      Real accum = 0;
      for (Int i = 0; i < 4; ++i) accum += cd[i]*q_ab_j[i];
      s_ab_j[d] = accum;
    }
    Real term2 = 0;
    for (Int d = 0; d < 3; ++d) term2 += s[d]*s_ab_j[d];
    term2 /= r2;
    for (Int d = 0; d < 3; ++d)
      s_ab_j[d] = (s_ab_j[d] - term2*s[d])/r;
  }
  for (Int d = 0; d < 3; ++d) s[d] /= r;
}

template <typename ConstVec3sT, typename Quad>
KOKKOS_INLINE_FUNCTION
void calc_sphere_to_ref (
  // The spherical quad containing the point.
  const ConstVec3sT& p, const Quad& e,
  // The point on the sphere.
  const Real q[3],
  // (a,b) in [-1,1]
  Real& a, Real& b,
  // Optional info output.
  Info* const info = nullptr,
  // Max number of iterations before returning with failure.
  const Int max_its = 10,
  // Tolerance for Newton iteration.
  const Real tol = 1e2*std::numeric_limits<Real>::epsilon())
{
  const Real tol2 = tol*tol;

  Real c[12];
  for (Int i = 0; i < 4; ++i)
    for (Int d = 0; d < 3; ++d)
      c[4*d + i] = p(e[i],d);

  Real rnorm2 = 1;
  a = b = 0;
  Int it = 0;
  for (Int it = 1; it < max_its; ++it) {
    Real s_ab[6], r[3], dab[2];
    calc_Jacobian(c, a, b, s_ab, r);
    for (Int d = 0; d < 3; ++d) r[d] -= q[d];
    rnorm2 = SphereGeometry::dot(r, r);
    if (rnorm2 <= tol2) break;
    solve_Jxr(s_ab, r, dab);
    a -= dab[0];
    b -= dab[1];
  }

  if (info) {
    info->success = rnorm2 <= tol2;
    info->n_iterations = it;
  }
}
} // namespace test

/* The jacobian when all nodes are moved. Use the interp formula for
   position. Let (a,b) in [-1,1]^2. Let va be the GLL basis values at a, and
   similarly for b. Let p(:,i) be the physical location of the i'th node. Then
       f(a,b) = sum_{i=1}^np sum_{j=1}^np p(:, np*(j-1) + i) va[i] vb[j]
       g(a,b) = norm(f(a,b))
       q = f(a,b) / g(a,b).
   Derivatives:
       q_a = f_a/g - (f g_a)/g^2
       g_a = g_f f_a
       g_f = 1/2 (f'f)^(-1/2) 2 f = f/norm(f)
   and similarly for q_b. In this case, we have
       f_a = sum_{i=1}^np sum_{j=1}^np p(:, np*(j-1) + i) va_a[i] vb[j],
   and va_a comes from the GLL polynomial.
*/
template <typename Cell>
static void calc_isoparametric_Jacobian (
  const AVec3s& ps, const Cell& cell, const Int np,
  const Real a, const Real b, Real J[6], Real r[3])
{
  for (int i = 0; i < 3; ++i) r[i] = 0;
  for (int i = 0; i < 6; ++i) J[i] = 0;
  {
    // Get the GLL values at ref coord (a,b).
    Real va[GLL::np_max], va_p[GLL::np_max];
    Real vb[GLL::np_max], vb_p[GLL::np_max];
    {
      GLL gll;
      gll.eval(np, a, va);
      gll.eval_derivative(np, a, va_p);
      gll.eval(np, b, vb);
      gll.eval_derivative(np, b, vb_p);
    }
    // Compute the isoparametric point in physical space corresponding to
    // (a,b), as well as derivatives of that function w.r.t. a and b.
    for (Int j = 0, k = 0; j < np; ++j)
      for (Int i = 0; i < np; ++i, ++k) {
        const auto& p = slice(ps, cell[k]);
        Real c = va[i]*vb[j];
        for (Int d = 0; d < 3; ++d) r[d]   += c*p[d];
        c = va_p[i]*vb[j];
        for (Int d = 0; d < 3; ++d) J[d]   += c*p[d];
        c = va[i]*vb_p[j];
        for (Int d = 0; d < 3; ++d) J[3+d] += c*p[d];
      }
  }
  // r'J
  Real rtJ[2] = {0};
  for (Int j = 0; j < 2; ++j) {
    const Real* const Jj = J + 3*j;
    for (Int i = 0; i < 3; ++i)
      rtJ[j] += r[i]*Jj[i];
  }
  const Real rnorm2 = siqk::SphereGeometry::norm2(r);
  const Real rnorm = std::sqrt(rnorm2);
  for (Int j = 0; j < 2; ++j) {
    Real* const Jj = J + 3*j;
    for (Int i = 0; i < 3; ++i)
      Jj[i] = (Jj[i] - r[i]*rtJ[j]/rnorm2)/rnorm;
  }
  for (Int i = 0; i < 3; ++i)
    r[i] /= rnorm;
}

template <typename Cell>
static Real calc_isoparametric_jacobian (
  const AVec3s& ps, const Cell& cell, const Int np,
  const Real a, const Real b)
{
  Real J[6], r[3];
  calc_isoparametric_Jacobian(ps, cell, np, a, b, J, r);
  SphereGeo::cross(J, J+3, r);
  return std::sqrt(SphereGeo::norm2(r));
}

template <typename ConstVec3sT, typename Cell>
void calc_isoparametric_sphere_to_ref (
  // The cell containing the point.
  const ConstVec3sT& p, const Cell& e, const Int np,
  // The point on the sphere.
  const Real q[3],
  // (a,b) in [-1,1]
  Real& a, Real& b,
  // Optional info output.
  siqk::sqr::Info* const info = nullptr,
  // Max number of iterations before returning with failure.
  const Int max_its = 10,
  // Tolerance for Newton iteration.
  const Real tol = 1e2*std::numeric_limits<Real>::epsilon())
{
  const Real tol2 = square(tol);
  Real rnorm2 = 1;
  a = b = 0;
  Int it = 0;
  for (it = 1; it <= max_its; ++it) { // Newton's method.
    Real J[6], r[3];
    calc_isoparametric_Jacobian(p, e, np, a, b, J, r);
    for (Int i = 0; i < 3; ++i) r[i] -= q[i];
    rnorm2 = siqk::SphereGeometry::norm2(r);
    if (rnorm2 <= tol2) break;
    Real dx[2];
    siqk::sqr::impl::solve_Jxr(J, r, dx);
    a -= dx[0];
    b -= dx[1];
  }
  if (info) {
    info->success = rnorm2 <= tol2;
    info->n_iterations = it;
  }
}

void Remapper::
init_isl_jacobian (const Mesh& m) {
  const Int dnn = nslices(m.dglln2cglln);
  resize(Jdep_, dnn);
  resize(Je_, dnn);
  const Int np = m.np, np2 = square(np);
  const Real* xnode;
  m.basis->get_x(m.np, xnode);
  for (Int i = 0; i < dnn; ++i) {
    const Int ti = i / np2;
    const auto& cell = slice(m.cgll_c2n, ti);
    const Int ni = i % np2;
    const Real ta = xnode[ni % np], tb = xnode[ni / np];
    if (md_) {
      // The correct, but increasingly expensive with np, method to compute
      // the density correction. Must use this if property preservation is on.
      Je_[i] = calc_isoparametric_jacobian(m.cgll_p, cell, m.np, ta, tb);
    } else {
      // Approximation to the density correction. When property preservation
      // is off, this is fine, as the tracer does not interact with the mass.
      // This approx means using just the corners, which means needing to
      // index advect_p properly.
      const Int corners[] = {cell[0], cell[np-1],
                             cell[np*np-1], cell[np*(np-1)]};
      Je_[i] = calc_jacobian(m.cgll_p, corners, ta, tb);
    }
  }
}

class FootprintTracker {
  Int ncell, np, np2;
  Array<Int> isl_nout, cisl_nout;
  Array<std::set<Int> > cisl_uniq_out;

public:
  typedef std::shared_ptr<FootprintTracker> Ptr;

  FootprintTracker (Int ncell_, Int np_)
    : ncell(ncell_), np(np_), np2(np*np),
      isl_nout(ncell, 0), cisl_nout(ncell, 0), cisl_uniq_out(ncell)
  {}

  void start () {
#   pragma omp parallel for
    for (Int i = 0; i < ncell; ++i) {
      isl_nout[i] = 0;
      cisl_uniq_out[i].clear();
    }
  }

  void record (Int tgt_gll, Int src_cell) {
    const Int ti = tgt_gll/np2;
    assert(ti >= 0 && ti < ncell);
    if (src_cell == ti) return;
#   pragma omp critical (FootprintTracker_record)
    {
      isl_nout[ti]++;
      cisl_uniq_out[ti].insert(src_cell);
    }
  }

  struct Stats {
    Int min, max, median;
    Real mean;
    void compute (Int* a, const Int n) {
      Int min = 1e3, max = -1;
      Real mean = 0;
#     pragma omp parallel
      {
#       pragma omp for reduction(min: min)
        for (Int i = 0; i < n; ++i) min = std::min(min, a[i]);
#       pragma omp for reduction(max: max)
        for (Int i = 0; i < n; ++i) max = std::max(max, a[i]);
#       pragma omp for reduction(+: mean)
        for (Int i = 0; i < n; ++i) mean += a[i];
      }
      mean /= n;
      std::nth_element(a, a + n/2 /* median or 1 higher */, a + n);
      median = a[n/2];
      this->min = min; this->max = max; this->mean = mean;
    }
  };

  void finish () {
#   pragma omp parallel for
    for (Int i = 0; i < ncell; ++i) {
      cisl_nout[i] = cisl_uniq_out[i].size();
      isl_nout[i] += 2*cisl_nout[i]; // min, max within source element
    }
    Stats isl;
    isl.compute(isl_nout.data(), ncell);
    printf("footprint> %2d %2d %4.1f %2d\n", isl.min, isl.median, isl.mean, isl.max);
  }
};

struct GlobalOnlyCedrData {
  std::vector<Real> Q_data, redistributed_mass, v;

  void set_min_cap (const size_t n) {
    if (Q_data.size() >= n) return;
    Q_data.resize(n); redistributed_mass.resize(n); v.resize(n);
  }
};

void glbl_only_pve (GlobalOnlyCedrData& d, const Int dnn, const Real* F_dgll,
                    const Real extra_mass, Real* rho_dgll, const bool mass_only = true) {
  if (mass_only) {
    Real vsum = 0;
    accumulate_threaded_bfb<1>([&] (const Int i, Real* a) { *a += F_dgll[i]; },
                               dnn, &vsum);
#   pragma omp parallel for
    for (Int i = 0; i < dnn; ++i) rho_dgll[i] += extra_mass/vsum;
    // For these test problems, we expect rho always to be >= 0, so we don't
    // impl a positivity filter here. But check that this is true.
#ifndef NDEBUG
    static bool rho_negative = false;
#   pragma omp parallel for
    for (Int i = 0; i < dnn; ++i)
      if ( ! rho_negative && rho_dgll[i] <= 0) {
        rho_negative = true;
        pr("WARNING rho <= 0:" pu(i) pu(rho_dgll[i]));
        assert(0);
      }
#endif
  }
}

void glbl_only_lcldyn (GlobalOnlyCedrData& d, const Int dnn,
                       const int ncell, const Real* qlim,
                       const Int* cn_src_cell, const Int* dgll2cglln,
                       const Real* F_dgll, const Real* rho_dgll,
                       Real& extra_mass, Real* Q_dgll, const bool mass_only = false) {
  if (mass_only) {
#   pragma omp parallel for
    for (Int i = 0; i < dnn; ++i)
      Q_dgll[i] += extra_mass/(4*M_PI);  
    return;
  }
  const Int np2 = dnn/ncell;
  Real adjust = 0;
  accumulate_threaded_bfb<1>(
    [&] (const Int i, Real* a) {
      const Int ci = cn_src_cell ? cn_src_cell[dgll2cglln[i]] : i/np2;
      const auto lo = qlim[2*ci]*rho_dgll[i], hi = qlim[2*ci+1]*rho_dgll[i];
      if (Q_dgll[i] < lo) {
        *a += F_dgll[i]*(Q_dgll[i] - lo);
        Q_dgll[i] = lo;
      } else if (Q_dgll[i] > hi) {
        *a += F_dgll[i]*(Q_dgll[i] - hi);
        Q_dgll[i] = hi;
      }
    },
    dnn, &adjust);
  extra_mass += adjust;
  if (extra_mass == 0) return;
  Real vsum = 0;
  if (extra_mass > 0) {
    accumulate_threaded_bfb<1>(
      [&] (const Int i, Real* a) {
        const Int ci = cn_src_cell ? cn_src_cell[dgll2cglln[i]] : i/np2;
        const auto hi = qlim[2*ci+1]*rho_dgll[i];
        if (Q_dgll[i] < hi) *a += F_dgll[i]*(hi - Q_dgll[i]);
      },
      dnn, &vsum);
  } else {
    accumulate_threaded_bfb<1>(
      [&] (const Int i, Real* a) {
        const Int ci = cn_src_cell ? cn_src_cell[dgll2cglln[i]] : i/np2;
        const auto lo = qlim[2*ci]*rho_dgll[i];
        if (Q_dgll[i] > lo) *a += F_dgll[i]*(Q_dgll[i] - lo);
      },
      dnn, &vsum);
  }
  const auto fac = extra_mass/vsum;
# pragma omp parallel for
  for (Int i = 0; i < dnn; ++i) {
    const Int ci = cn_src_cell ? cn_src_cell[dgll2cglln[i]] : i/np2;
    const auto lo = qlim[2*ci]*rho_dgll[i], hi = qlim[2*ci+1]*rho_dgll[i];
    const Real v = ((extra_mass > 0) ?
                    ((Q_dgll[i] < hi) ? (hi - Q_dgll[i]) : 0) :
                    ((Q_dgll[i] > lo) ? (Q_dgll[i] - lo) : 0));
    Q_dgll[i] += fac*v;
  }
}

// Adjust non-boundary GLL nodes so that r_t has the same mass as r_s.
static void
adjust_mass_nonbdy (const Real* w_s, const Real* r_s, const Int np2_s,
                    const Real* w_t, Real* r_t, const Int np_t) {
  const Int np2_t = square(np_t);
  Real m_s = 0, m_t = 0;
  for (Int i = 0; i < np2_s; ++i) m_s += w_s[i]*r_s[i];
  for (Int i = 0; i < np2_t; ++i) m_t += w_t[i]*r_t[i];
  const Real delta = m_s - m_t;
  Real den = 0;
  for (Int i = 1; i < np_t-1; ++i)
    for (Int j = 1; j < np_t-1; ++j)
      den += w_t[i*np_t + j];
  const Real dr = delta/den;
  for (Int i = 1; i < np_t-1; ++i)
    for (Int j = 1; j < np_t-1; ++j)
      r_t[i*np_t + j] += dr;
#ifndef NDEBUG
  for (Int i = 0; i < np2_t; ++i) assert(r_t[i] >= 0);
#endif
}

struct TendencyData {
  typedef std::shared_ptr<TendencyData> Ptr;
  static constexpr Real lcl_cdr_relax_fac = 1e-2;
  bool dss_source; // v-mesh DSS is an apparent source
  bool cdr_on_src_tends;
  const Real* dq; // mixing ratio
  Int idx_beg, idx_end; // tracers [idx_beg,idx_end) have tendencies
  Array<Real> q_min, q_max;

  TendencyData () : cdr_on_src_tends(false) {}
};

class Remapper::IslImpl {
  PRefineData::Ptr p;
  AVec3s p_fine;
  bool first_call;
  MeshInterpolator::Ptr mesh_interp;
  Transferer2D::Ptr t2d_c2f, t2d_f2c, t2d_c2f_rho;
  Array<Real> ls_f2c_op;
  
public:  
  Int idxs[2];
  Array<Real> tracer[2], rho[2];
  std::shared_ptr<spf::MassRedistributor> mrd_c, mrd_f; // (c)oarse and (f)ine
  FootprintTracker::Ptr footprint;
  D2Cer::Ptr d2cer_f;
  GlobalOnlyCedrData gocd;
  mutable Array<Real> wrk;
  TendencyData::Ptr tend;

  IslImpl (const PRefineData::Ptr& pr_)
    : p(pr_), first_call(true)
  {
    idxs[0] = 0; idxs[1] = 1;
    assert(p->m_t->basis);
    tend = std::make_shared<TendencyData>();
    tend->cdr_on_src_tends = p->md_t && ! Filter::is_nodal(p->md_t->global_type());
    if (p->experiment > 0) {
      assert(p->m_v->basis && p->m_f->basis);
      const auto gll_basis = Basis::create(Basis::Type::gll);
      SIQK_THROW_IF( ! p->m_v->basis->gll_nodes(),
                    "The coarse mesh should always have GLL nodes.");
      // Move the mesh using GLL basis.
      mesh_interp = std::make_shared<MeshInterpolator>(
        p->m_v->np, gll_basis, p->m_f->np, p->m_f->basis);
      // Interpolate from coarse to fine using GLL basis.
      t2d_c2f = std::make_shared<Interpolator2D>(
        p->m_v->np, gll_basis, p->m_f->np, p->m_f->basis);
      t2d_f2c = std::make_shared<Interpolator2D>(
        p->m_f->np, p->m_f->basis, p->m_v->np, gll_basis);
      t2d_c2f_rho = t2d_c2f;
    } else {
      SIQK_THROW_IF( ! p->m_t->basis->gll_nodes(),
                    "For prefine 0, the nodes must be GLL points.");
    }
  }

  void track_footprint (const bool track) {
    footprint = std::make_shared<FootprintTracker>(
      nslices(p->m_t->geo_c2n), p->m_f->np);
  }

  const AVec3s& interp_p (const AVec3s& p_coarse) {
    if (nslices(p_fine) == 0)
      p_fine = AVec3s(nslices(p->m_f->cgll_p));
    mesh_interp->apply(*p->m_v, *p->m_f, p_coarse, p_fine);
    return p_fine;
  }

  void init_tracer_arrays (const Int ntracers) {
    if (p->experiment == 0) return;
    for (Int i = 0; i < 2; ++i)
      tracer[i].optclear_and_resize(ntracers*len(p->m_f->dglln2cglln));
    for (Int i = 0; i < 2; ++i)
      rho[i].optclear_and_resize(ntracers*len(p->m_f->dglln2cglln));
  }

  bool start_time_step (const Int ntracers, const Real* srcterm_tendency,
                        const Int srcterm_beg, const Int srcterm_end,
                        const pg::PhysgridData::Ptr& pg) {
    assert(first_call || tend->dq == srcterm_tendency);
    if ( ! first_call) return false;
    init_tracer_arrays(ntracers);
    {
      tend->dq = srcterm_tendency;
      tend->idx_beg = srcterm_beg;
      tend->idx_end = srcterm_end;
      tend->dss_source = ! pg && tend->idx_beg < ntracers;
      if (srcterm_tendency || tend->dss_source) {
        const auto ncell = nslices(p->m_t->geo_c2n);
        tend->q_min.optclear_and_resize(ntracers*ncell);
        tend->q_max.optclear_and_resize(ntracers*ncell);
      }
    }
    first_call = false;
    return true;
  }

  void interp_rho (const Real* const src_density, Real* const rho,
                   const bool continuous) {
    const Int ne = nslices(p->m_v->geo_c2n);
    const Int np_src = p->m_v->np, np_tgt = p->m_f->np;
    const Int np2_c = square(np_src);
    const Int np2_f = square(np_tgt);
    const bool homme_mass = Dmc::use_homme_mass(p->rd_f->dmc());
    const Real* w_src = p->rd_v->dgbfi_mass().data();
    const Real* w_tgt = p->rd_f->dgbfi_mass().data();
    const Real* jac_src = p->rd_v->Jt().data();
    const Real* jac_tgt = p->rd_f->Jt().data();
#   pragma omp parallel for
    for (Int ie = 0; ie < ne; ++ie) {
      if (homme_mass) {
        const Real* src_ie = src_density + ie*np2_c;
        Real* tgt_ie = rho + ie*np2_f;
        if ( ! continuous) {
          Real src[GLL::np_max*GLL::np_max];
          const Real* jac_src_ie = jac_src + ie*np2_c;
          for (Int i = 0; i < np2_c; ++i) src[i] = src_ie[i]*jac_src_ie[i];
          t2d_c2f->apply(src, tgt_ie);
          const Real* jac_tgt_ie = jac_tgt + ie*np2_f;
          for (Int i = 0; i < np2_f; ++i) tgt_ie[i] /= jac_tgt_ie[i];
        } else {
          t2d_c2f->apply(src_ie, tgt_ie);
          const Real* w_src_ie = w_src + ie*np2_c;
          const Real* w_tgt_ie = w_tgt + ie*np2_f;
          adjust_mass_nonbdy(w_src_ie, src_ie, np2_c, w_tgt_ie, tgt_ie, np_tgt);
        }
#if 0
        {
          const Real* w_src = p->rd_v->dgbfi_mass().data() + ie*np2_c;
          const Real* w_tgt = p->rd_f->dgbfi_mass().data() + ie*np2_f;
          Real m_src = 0, m_tgt = 0;
          for (Int i = 0; i < np2_c; ++i) m_src += w_src[i]*src_ie[i];
          for (Int i = 0; i < np2_f; ++i) m_tgt += w_tgt[i]*tgt_ie[i];
          if (reldif(m_src,m_tgt) > 2e-15)
            pr(puf(ie) pu(m_src) pu(m_tgt) pu(reldif(m_src,m_tgt)));
        }
#endif
      } else {
        t2d_c2f->apply(src_density + ie*np2_c, rho + ie*np2_f);
      }
    }
  }

  void interp_q_to_t_mesh (const Real* const tracer_from, Real* const tracer_to,
                           const Int ntracers) {
    const Int ne = nslices(p->m_v->geo_c2n);
    const Int np_v2 = square(p->m_v->np), np_t2 = square(p->m_f->np);
    assert(np_v2 == square(t2d_c2f->get_np_from()));
    assert(np_t2 == square(t2d_c2f->get_np_to()));
#   pragma omp parallel
    {
      for (Int ti = 0; ti < ntracers; ++ti) {
#       pragma omp for
        for (Int ie = 0; ie < ne; ++ie)
          t2d_c2f->apply(tracer_from + (ti*ne + ie)*np_v2,
                         tracer_to + (ti*ne + ie)*np_t2);
      }
    }    
  }
  
  void interp_q_to_v_mesh (const Real* const tracer_from, Real* const tracer_to,
                           const Int ntracers) {
    const Int ne = nslices(p->m_v->geo_c2n);
    const Int np_v2 = square(p->m_v->np), np_t2 = square(p->m_f->np);
    assert(np_t2 == square(t2d_f2c->get_np_from()));
    assert(np_v2 == square(t2d_f2c->get_np_to()));
#   pragma omp parallel
    {
      for (Int ti = 0; ti < ntracers; ++ti) {
#       pragma omp for
        for (Int ie = 0; ie < ne; ++ie)
          t2d_f2c->apply(tracer_from + (ti*ne + ie)*np_t2,
                         tracer_to + (ti*ne + ie)*np_v2);
      }
    }
  }

  void project_q_to_v_mesh (const Real* const rho_from, const Real* const rho_to,
                            const Real* const q_from, Real* const q_to,
                            const Int ntracers) {
    const Int ncell = nslices(p->m_v->geo_c2n);
    const Int np_v = p->m_v->np, np2_v = square(np_v), np_t = p->m_f->np,
      np2_t = square(np_t);
    if (ls_f2c_op.empty()) {
      Array<Real> A(np2_v*np2_v), op(np2_v*np2_t);
      Basis::compute_mass_matrix_2d(*p->m_v->basis, np_v, *p->m_v->basis, np_v,
                                    A.data());
      Basis::compute_mass_matrix_2d(*p->m_t->basis, np_t, *p->m_v->basis, np_v,
                                    op.data());
      const Int e = form_ls_op(np2_v, np2_t, A.data(), op.data());
      SIQK_THROW_IF(e != 0, "Couldn't form project_q_to_v_mesh op");
      ls_f2c_op.optclear_and_resize(np2_v*np2_t);
      transpose(np2_v, np2_t, op.data(), ls_f2c_op.data());
    }
#   pragma omp parallel
    {
      const Real* const jacs_t = p->rd_f->Jt().data();
      const Real* const jacs_v = p->rd_v->Jt().data();
      for (Int ti = 0; ti < ntracers; ++ti) {
#       pragma omp for
        for (Int ie = 0; ie < ncell; ++ie) {
          const Real* jac_t = jacs_t + ie*np2_t;
          const Real* rho_t = rho_from + ie*np2_t;
          const Real* q_t = q_from + (ti*ncell + ie)*np2_t;
          const Real* jac_v = jacs_v + ie*np2_v;
          const Real* rho_v = rho_to + ie*np2_v;
          Real* q_v = q_to + (ti*ncell + ie)*np2_v;
          Real src[GLL::np_max*GLL::np_max], tgt[GLL::np_max*GLL::np_max];
          for (Int i = 0; i < np2_t; ++i) src[i] = jac_t[i]*rho_t[i]*q_t[i];
          matmult_rcc(np2_v, 1, np2_t, 1, ls_f2c_op.data(), src, 0, tgt);
          for (Int i = 0; i < np2_v; ++i) q_v[i] = tgt[i]/(jac_v[i]*rho_v[i]);
#if 0
          {
            const Real* w_t = p->rd_f->dgbfi_mass().data() + ie*np2_t;
            const Real* w_v = p->rd_v->dgbfi_mass().data() + ie*np2_v;
            Real m_t = 0, m_v = 0;
            for (Int i = 0; i < np2_t; ++i) m_t += w_t[i]*rho_t[i]*q_t[i];
            for (Int i = 0; i < np2_v; ++i) m_v += w_v[i]*rho_v[i]*q_v[i];
            if (reldif(m_t,m_v) > 1e-13)
              pr(puf(ti) pu(ie) pu(m_t) pu(m_v) pu(reldif(m_t,m_v)));
          }
#endif
        }
      }
    }
  }

  void limit_tracers (
    const Mesh& m_from, const Mesh& m_to,
    const RemapData& rd_from, const RemapData& rd_to,
    const MonoData& md_to,
    const Real* const rho_from, const Real* const tracer_from,
    const Real* const rho_to, Real* const tracer_to, const Int ntracers)
  {
    const Int np_to = m_to.np, np2_to = square(np_to);
    const Int np_from = m_from.np, np2_from = square(np_from);
    const auto dgbfi_mass_to = rd_to.dgbfi_mass().data();
    const auto dgbfi_mass_from = rd_from.dgbfi_mass().data();
    const Int ncell = nslices(m_to.geo_c2n);
    for (Int ti = 0; ti < ntracers; ++ti) {
#     pragma omp parallel for
      for (Int ci = 0; ci < ncell; ++ci) {
        const Real* const F_from = dgbfi_mass_from + ci*np2_from;
        const Real* const r_from = rho_from + ci*np2_from;
        const Real* const t_from = tracer_from + (ti*ncell + ci)*np2_from;
        const Real* const F_to = dgbfi_mass_to + ci*np2_to;
        const Real* const r_to = rho_to + ci*np2_to;
        Real* const t_to = tracer_to + (ti*ncell + ci)*np2_to;
        Real Qm_to = 0, Qm_from = 0, q_min = t_from[0], q_max = q_min;
        for (Int i = 0; i < np2_to; ++i) {
          t_to[i] *= r_to[i];
          Qm_to += F_to[i]*t_to[i];
        }
        for (Int i = 0; i < np2_from; ++i)
          Qm_from += F_from[i]*r_from[i]*t_from[i];
        for (Int i = 1; i < np2_from; ++i)
          q_min = std::min(q_min, t_from[i]);
        for (Int i = 1; i < np2_from; ++i)
          q_max = std::max(q_max, t_from[i]);
        // This function is called only on feasible problems.
        md_to.limit_tracer(ci, q_min, q_max, r_to, t_to, false /* expand */,
                           Qm_from - Qm_to, true /* nocheck */);
        for (Int i = 0; i < np2_to; ++i)
          t_to[i] /= r_to[i];
      }
    }
  }

  void transfer_q_to_t_mesh (const bool run_cdr, const Int ntracers,
                             const Real* const rho_from, const Real* const tracer_from,
                             const Real* const rho_to, Real* const tracer_to) {
    interp_q_to_t_mesh(tracer_from, tracer_to, ntracers);
    if (run_cdr) {
      // Limit q on v mesh. This is local only.
      limit_tracers(*p->m_v, *p->m_f, *p->rd_v, *p->rd_f, *p->md_f,
                    rho_from, tracer_from,
                    rho_to, tracer_to, ntracers);
    }
  }

  void transfer_q_to_v_mesh (const bool run_cdr, const Int ntracers,
                             const Real* const rho_from, const Real* const tracer_from,
                             const Real* const rho_to, Real* const tracer_to) {
    // Projecting Q is worse for toy chem diagnostic measured on v grid.
    interp_q_to_v_mesh(tracer_from, tracer_to, ntracers);
    //project_q_to_v_mesh(rho_from, rho_to, tracer_from, tracer_to, ntracers);
    if (run_cdr) {
      // Limit q on v mesh. This is local only.
      limit_tracers(*p->m_f, *p->m_v, *p->rd_f, *p->rd_v, *p->md_v,
                    rho_from, tracer_from,
                    rho_to, tracer_to, ntracers);
    }
  }

  void transfer_to_physgrid (pg::PhysgridData& pg, const bool run_cdr,
                             const Int ntracers, const Real* const rho_t,
                             const Real* const q_t) {
    const auto md = *pg.ops->mesh_data;
    const auto op = *pg.ops->gll2fv;
    const Int np2 = square(pg.ops->np), nf2 = square(pg.ops->nphys), ncell = md.nelem;
#   pragma omp parallel for
    for (Int ci = 0; ci < ncell; ++ci)
      for (Int ti = 0; ti < ntracers; ++ti) {
        op.remap(&md.gll_metdet[ci*np2], &md.fv_metdet[ci*nf2], pg.limiter,
                 &rho_t[ci*np2], &q_t[ncell*np2*ti + ci*np2],
                 &pg.rho[ci*nf2], &pg.q[ncell*nf2*ti + ci*nf2],
                 ti == 0 /* remap rho only the first time */);
#if 0
        {
          const Real wt_pg = 1.0/nf2;
          const Real* jac_pg = &md.fv_metdet[ci*nf2];
          const Real* wtjac_gll = p->rd_f->dgbfi_mass().data() + ci*np2;
          const Real* rho_pg = &pg.rho[ci*nf2];
          const Real* rho_gll = &rho_t[ci*np2];
          const Real* q_pg = &pg.q[ncell*nf2*ti + ci*nf2];
          const Real* q_gll = &q_t[ncell*np2*ti + ci*np2];
          Real Qm_pg = 0, Qm_gll = 0;
          for (Int i = 0; i < nf2; ++i) Qm_pg += wt_pg*jac_pg[i]*rho_pg[i]*q_pg[i];
          for (Int i = 0; i < np2; ++i) Qm_gll += wtjac_gll[i]*rho_gll[i]*q_gll[i];
          const Real re = reldif(Qm_pg, Qm_gll);
          if (re > 1e-13) pr("transfer" pu(ti) pu(ci) pu(Qm_pg) pu(Qm_gll) pu(re));
        }
#endif
      }
  }

  void add_tendencies_from_vgrid (
    const Real* const rho_from, const Real* const tracer_from,
    const Real* const rho_to, Real* const tracer_to,
    const bool run_cdr, const Int ntracers, TendencyData& td, const Real* dq)
  {
    if (td.idx_beg >= ntracers) return;
    const auto& m_from = *p->m_v;
    const auto& m_to = *p->m_f;
    const auto& rd_to = *p->rd_f;
    const auto md_to = p->md_f;
    const Int np_from = m_from.np, np2_from = square(np_from);
    const Int np_to = m_to.np, np2_to = square(np_to);
    const auto jacobian_from = p->rd_v->Jt();
    const auto jacobian_to   = p->rd_f->Jt();
    const Int ncell = nslices(m_to.geo_c2n);
    const bool homme_mass = Dmc::use_homme_mass(rd_to.dmc());
    const Int dnn_from = ncell*np2_from;
    if (td.dss_source && wrk.size() < (size_t) dnn_from)
      wrk.optclear_and_resize(dnn_from);
    for (Int ti = 0; ti < ntracers; ++ti) {
      const bool have_tendency = (dq && ti >= td.idx_beg && ti < td.idx_end);
      if ( ! (td.dss_source || have_tendency)) continue;
      if (td.dss_source) {
        // q on the v mesh prior to the DSS.
        transfer_q_to_v_mesh(run_cdr, 1, rho_to, tracer_to + np2_to*ncell*ti,
                             rho_from, wrk.data());
      }
#     pragma omp parallel for
      for (Int ci = 0; ci < ncell; ++ci) {
        const Real* const r_from = rho_from + ci*np2_from;
        const Real* const J_from = jacobian_from.data() + ci*np2_from;
        const Real* const q_from = tracer_from + (ti*ncell + ci)*np2_from;
        const Real* const r_to = rho_to + ci*np2_to;
        const Real* const J_to = jacobian_to.data() + ci*np2_to;
        Real* const q_to = tracer_to + (ti*ncell + ci)*np2_to;
        const Real* const dq_from = have_tendency ?
          (dq + ((ti - td.idx_beg)*ncell + ci)*np2_from) :
          nullptr;
        Real dQ_from[GLL::np_max*GLL::np_max] = {0};
        if (td.dss_source) {
          // Account for apparent source due to DSS on v mesh.
          const Real* const q_v = wrk.data() + ci*np2_from;
          for (Int i = 0; i < np2_from; ++i)
            dQ_from[i] = r_from[i]*(q_from[i] - q_v[i]);
        }
        if (have_tendency) // actual source term
          for (Int i = 0; i < np2_from; ++i)
            dQ_from[i] += r_from[i]*dq_from[i];
        // Interpolate density, which is conservative.
        Real Q_to[GLL::np_max*GLL::np_max];
        if (homme_mass) {
          Real src[GLL::np_max*GLL::np_max];
          for (Int i = 0; i < np2_from; ++i) src[i] = dQ_from[i]*J_from[i];
          t2d_c2f->apply(src, Q_to);
          for (Int i = 0; i < np2_to; ++i) Q_to[i] /= J_to[i]; 
        } else {
          t2d_c2f->apply(dQ_from, Q_to);
        }
        // Get bounds.
        Real q_min = -1, q_max = -1;
        for (Int i = 0; i < np2_from; ++i) {
          const Real q = have_tendency ? (q_from[i] + dq_from[i]) : q_from[i];
          if (i == 0)
            q_min = q_max = q;
          else {
            q_min = std::min(q_min, q);
            q_max = std::max(q_max, q);
          }
        }
        // Augment with current fine-mesh bounds to assure that if dQ = 0, then
        // the fine-mesh tracer is not modified.
        for (Int i = 0; i < np2_to; ++i) q_min = std::min(q_min, q_to[i]);
        for (Int i = 0; i < np2_to; ++i) q_max = std::max(q_max, q_to[i]);
        if ( ! td.cdr_on_src_tends) {
          const Int os = ti*ncell;
          td.q_min[os + ci] = q_min;
          td.q_max[os + ci] = q_max;
        }
        // Update fine-mesh tracer.
        if (run_cdr) {
          for (Int i = 0; i < np2_to; ++i) Q_to[i] += r_to[i]*q_to[i];
          const Real del = td.cdr_on_src_tends ? 0 :
            TendencyData::lcl_cdr_relax_fac*(q_max - q_min);
          Int info =
            md_to->limit_tracer(ci, q_min - del, q_max + del, r_to, Q_to,
                                ! td.cdr_on_src_tends /* expand */, 0,
                                true /* nocheck */ , md_to->local_type());
          for (Int i = 0; i < np2_to; ++i) q_to[i] = Q_to[i]/r_to[i];
          if (false && info < 0) {
            Real new_max = 0, new_min = 1;
            for (Int i = 0; i < np2_to; ++i) new_max = std::max(new_max, q_to[i]);
            for (Int i = 0; i < np2_to; ++i) new_min = std::min(new_min, q_to[i]);
            pr("tend ADD" pu(ti) pu(ci) pu(info));
            if (new_min < q_min) pr(puf(q_min) pu(new_min) pu(q_min-new_min));
            if (new_max > q_max) pr(puf(q_max) pu(new_max) pu(new_max-q_max));
          }
        } else {
          for (Int i = 0; i < np2_to; ++i) q_to[i] += Q_to[i]/r_to[i];
        }
      }
    }
  }

  void add_tendencies_from_physgrid (
    const Real* const rho, Real* const tracer,
    const bool run_cdr, const Int ntracers, TendencyData& td, const Real* dq,
    const pg::PhysgridData& pg, D2Cer::Ptr d2cer = nullptr)
  {
    assert(Dmc::use_homme_mass(p->rd_f->dmc()));
    // td.dss_source is irrelevant here because there is no DSS when the tracers
    // are on the physgrid.
    if ( ! dq) return;
    const Int nphys2 = square(pg.ops->nphys);
    const Int ncell = pg.ops->mesh_data->nelem;
    const Int np = pg.ops->np, np2 = square(np);
    const auto& md = *pg.ops->mesh_data;
    const auto& op = *pg.ops->fv2gll;
    const auto md_to = p->md_f;
    if (wrk.size() < size_t(2*ncell)) wrk.optclear_and_resize(2*ncell);
    for (Int ti = td.idx_beg; ti < td.idx_end; ++ti)
#   pragma omp parallel
    {
#     pragma omp for
      for (Int ci = 0; ci < ncell; ++ci) {
        const Real* const q_pg = pg.q.data() + (ti*ncell + ci)*nphys2;
        const Real* const dq_pg = dq + ((ti - td.idx_beg)*ncell + ci)*nphys2;
        Real q_min = -1, q_max = -1;
        for (Int i = 0; i < nphys2; ++i) {
          const Real q = q_pg[i] + dq_pg[i];
          if (i == 0)
            q_min = q_max = q;
          else {
            q_min = std::min(q_min, q);
            q_max = std::max(q_max, q);
          }
        }
        wrk[2*ci  ] = q_min;
        wrk[2*ci+1] = q_max;
      }
      const Int os = ti*ncell;
#     pragma omp for
      for (Int ci = 0; ci < ncell; ++ci) {
        Real q_min = wrk[2*ci], q_max = wrk[2*ci+1];
        for (Int j = md.geo_c2cnbrs_ptr[ci]; j < md.geo_c2cnbrs_ptr[ci+1]; ++j) {
          q_min = std::min(q_min, wrk[2*md.geo_c2cnbrs[j]  ]);
          q_max = std::max(q_max, wrk[2*md.geo_c2cnbrs[j]+1]);
        }
        // Augment with current fine-mesh bounds to assure that if dQ = 0, then
        // the fine-mesh tracer is not modified.
        const Real* const q_gll = tracer + (ncell*ti + ci)*np2;
        for (Int i = 0; i < np2; ++i) q_min = std::min(q_min, q_gll[i]);
        for (Int i = 0; i < np2; ++i) q_max = std::max(q_max, q_gll[i]);
        td.q_min[os + ci] = q_min; td.q_max[os + ci] = q_max;
      }
#     pragma omp for
      for (Int ci = 0; ci < ncell; ++ci) {
        const Real* const dq_pg = dq + ((ti - td.idx_beg)*ncell + ci)*nphys2;
        Real dq_gll[GLL::np_max*GLL::np_max];
        op.remap(&md.gll_metdet[ci*np2], &md.fv_metdet[ci*nphys2],
                 Limiter::none, -1, -1,
                 &pg.rho[ci*nphys2], dq_pg,
                 const_cast<Real*>(&rho[ci*np2]), dq_gll,
                 false /* don't remap rho */);
        Real* const q_gll = tracer + (ncell*ti + ci)*np2;
        for (Int i = 0; i < np2; ++i) q_gll[i] += dq_gll[i];
        if (run_cdr) {// && td.cdr_on_src_tends) {
          const Real* rho_gll = &rho[ci*np2];
          for (Int i = 0; i < np2; ++i) q_gll[i] *= rho_gll[i];
          // Relax the constraints and expand them, if necessary, when we're
          // not solving local problems strictly.
          const Real del = td.cdr_on_src_tends ? 0 :
            TendencyData::lcl_cdr_relax_fac*(td.q_max[os + ci] - td.q_min[os + ci]);
          md_to->limit_tracer(ci, td.q_min[os + ci] - del, td.q_max[os + ci] + del,
                              &rho[ci*np2], q_gll, ! td.cdr_on_src_tends, 0, true);
          for (Int i = 0; i < np2; ++i) q_gll[i] /= rho_gll[i];
        }
#if 0
        {
          const Real wt_pg = 1.0/nphys2;
          const Real* jac_pg = &md.fv_metdet[ci*nphys2];
          const Real* wtjac_gll = p->rd_f->dgbfi_mass().data() + ci*np2;
          const Real* rho_pg = &pg.rho[ci*nphys2];
          const Real* q_pg = &pg.q[ncell*nphys2*ti + ci*nphys2];
          const Real* rho_gll = &rho[ci*np2];
          Real Qm_pg = 0, Qm_gll = 0;
          for (Int i = 0; i < nphys2; ++i) Qm_pg += wt_pg*jac_pg[i]*rho_pg[i]*(q_pg[i] + dq_pg[i]);
          for (Int i = 0; i < np2; ++i) Qm_gll += wtjac_gll[i]*rho_gll[i]*q_gll[i];
          const Real re = reldif(Qm_pg, Qm_gll);
          if (re > 1e-11) pr("tend" pu(ti) pu(ci) pu(Qm_pg) pu(Qm_gll) pu(re));
          if (md_to)
            for (Int i = 0; i < np2; ++i)
              assert(md_to->dgbfi_mass()[ci*np2 + i] == wtjac_gll[i]);
        }
#endif
      }
    }
  }

  void add_tendencies (
    const Real* const rho_from, const Real* const tracer_from,
    const Real* const rho_to, Real* const tracer_to,
    const bool run_cdr, const Int ntracers, TendencyData& td, const Real* dq,
    const pg::PhysgridData::Ptr& pg)
  {
    if (pg)
      add_tendencies_from_physgrid(rho_to, tracer_to,
                                   run_cdr, ntracers, td, dq, *pg);
    else
      add_tendencies_from_vgrid(rho_from, tracer_from, rho_to, tracer_to,
                                run_cdr, ntracers, td, dq);
  }

  void dss_q (const Real* const rho, Real* const q, const Int ntracers) {
    if (! d2cer_f)
      d2cer_f = std::make_shared<D2Cer>(p->m_f->dglln2cglln, p->rd_f->dgbfi_mass());
    const Int dnn = d2cer_f->get_dnn();
    if (wrk.size() < size_t(dnn)) wrk.optclear_and_resize(dnn);
    for (Int ti = 0; ti < ntracers; ++ti) {
      Real* t = q + ti*dnn;
      d2cer_f->dss_q(rho, t, wrk.data());
    }
  }

  void dss_q (const Real* const rho, Real* const tracer, const Int tb, const Int te) {
    if ( ! d2cer_f)
      d2cer_f = std::make_shared<D2Cer>(p->m_f->dglln2cglln, p->rd_f->dgbfi_mass());
    const Int dnn = d2cer_f->get_dnn();
    if (wrk.size() < size_t(dnn)) wrk.optclear_and_resize(dnn);
    for (Int ti = tb; ti < te; ++ti) {
      Real* t = tracer + ti*dnn;
      d2cer_f->dss_q(rho, t, wrk.data());
    }
  }

  void finish_time_step () {
    std::swap(idxs[0], idxs[1]);
  }
};

void Remapper::track_footprint (const bool track) {
  isl_impl_->track_footprint(track);
}

void Remapper::make_c2d_relations (
  const Int cnn, const AIdxArray& dglln2cglln,
  C2DRelations& c2d)
{
  const Int dnn = nslices(dglln2cglln);
  resize(c2d.ptr, cnn+1);
  resize(c2d.r, dnn);
  copy(c2d.ptr, 0);
  // Count.
  for (Int i = 0; i < dnn; ++i) {
    const Int ci = dglln2cglln[i];
    if (ci == cnn-1) continue;
    ++c2d.ptr[ci+2];
  }
  // Cumsum.
  for (Int i = 3; i <= cnn; ++i)
    c2d.ptr[i] += c2d.ptr[i-1];
  // Fill.
  for (Int i = 0; i < dnn; ++i)
    c2d.r[ c2d.ptr[ dglln2cglln[i]+1 ]++ ] = i;
#if 1
  // Test.
  assert(c2d.ptr[0] == 0);
  assert(c2d.ptr[cnn] == dnn);
  for (Int i = 0; i < dnn; ++i) {
    const Int ci = dglln2cglln[i];
    bool found = false;
    for (Int j = c2d.ptr[ci]; j < c2d.ptr[ci+1]; ++j)
      if (c2d.r[j] == i) {
        found = true;
        break;
      }
    assert(found);
  }
#endif
}

void Remapper::
init_isl () {
  isl_impl_ = std::make_shared<IslImpl>(pr_);
  const auto& m = *(pr_->experiment == 0 ? pr_->m_f : pr_->m_v);
  const Int cnn = nslices(m.cgll_p);
  make_c2d_relations(nslices(m.cgll_p), m.dglln2cglln, c2d_v_);
  Int cnn_f = 0;
  if (pr_->experiment > 0) {
    cnn_f = nslices(pr_->m_f->cgll_p);
    make_c2d_relations(cnn_f, pr_->m_f->dglln2cglln, c2d_f_);
  }
  resize(cn_src_cell_, std::max(cnn, cnn_f));
  // independent of mesh
  resize(q_data_, ncell_, 2);
  init_isl_jacobian(m);
}

// Always d2c rho because the density factor is discontinuous at the element
// edges.
static void dss_rho (D2Cer& d2cer, Real* rho, Array<Real>& wrk) {
  const auto dnn = d2cer.get_dnn();
  if (wrk.size() < (size_t) dnn) wrk.optclear_and_resize(dnn);
  d2cer.dss(rho, wrk.data());
}

static Int
find_src_cell (const Mesh& m, const RemapData& rd, const Int ne,
               const Int cni /* target node, continuous ID */,
               const AVec3s& advected_p) {
  // Solve for departure point.
  const auto p = slice(advected_p, cni);
  // Normalize p for much faster calc_sphere_to_ref.
  Real p_sphere[3];
  for (Int d = 0; d < 3; ++d) p_sphere[d] = p[d];
  siqk::SphereGeometry::normalize(p_sphere);
  Int ci;
  if (m.nonuni) {
    SearchFunctor sf(m, p_sphere);
    rd.octree().apply(sf.bb, sf); // same regardless of mesh
    ci = sf.ci_hit;
  } else {
    const auto& gr = m.grid_rotation;
    ci = mesh::get_cell_idx(ne, gr.angle, gr.R, p_sphere[0], p_sphere[1],
                            p_sphere[2]);
  }
  return ci;
}

static void
calc_departure_data (const Mesh& m, const RemapData& rd, const Int ne,
                     const Int cni /* target node, continuous ID */,
                     const AVec3s& advected_p,
                     Int& ci, bool calc_ci, Real* va, Real* vb, Real& a, Real& b) {
  // Solve for departure point.
  const auto p = slice(advected_p, cni);
  // Normalize p for much faster calc_sphere_to_ref.
  Real p_sphere[3];
  for (Int d = 0; d < 3; ++d) p_sphere[d] = p[d];
  siqk::SphereGeometry::normalize(p_sphere);
  if (calc_ci) {
    if (m.nonuni) {
      SearchFunctor sf(m, p_sphere);
      rd.octree().apply(sf.bb, sf); // same regardless of mesh
      ci = sf.ci_hit;
    } else {
      const auto& gr = m.grid_rotation;
      ci = mesh::get_cell_idx(ne, gr.angle, gr.R, p_sphere[0], p_sphere[1],
                              p_sphere[2]);
    }
  }
  assert(ci >= 0 && ci < nslices(m.dglln2cglln)/square(m.np));
  const auto& cell = slice(m.geo_c2n, ci);
  siqk::sqr::Info info;
  sphere2ref::calc_sphere_to_ref(m.geo_p, cell, p_sphere, a, b, &info);
  // Eval basis functions.
  m.basis->eval(m.np, a, va);
  m.basis->eval(m.np, b, vb);
}

static void
calc_jacobian_departure (const Mesh& m, const bool isoparametric,
                         const AVec3s& advected_p, Real* Jdep) {
  const Int dnn = nslices(m.dglln2cglln);
  const Int np = m.np, np2 = square(np);
  const Real* xnode;
  m.basis->get_x(m.np, xnode);
  ompparfor for (Int i = 0; i < dnn; ++i) {
    const Int ti = i / np2;
    const auto& cell = slice(m.cgll_c2n, ti);
    const Int ni = i % np2;
    const Real ta = xnode[ni % np], tb = xnode[ni / np];
    Real Jd;
    if (isoparametric) {
      // This is the case of interest and runs when property preservation is
      // on. For the case of p-refinement, it runs with np=npv=4, not npt.
      Jd = calc_isoparametric_jacobian(advected_p, cell, np, ta, tb);
    } else {
      // This happens only when property preservation is off, in which case
      // density doesn't couple to the mixing ratios.
      const Int corners[] = {cell[0], cell[np-1],
                             cell[np*np-1], cell[np*(np-1)]};
      Jd = calc_jacobian(advected_p, corners, ta, tb);
    }
    Jdep[i] = Jd;
  }
}

void Remapper::
interp (const Mesh& m, const C2DRelations& c2d, const AVec3s& advected_p,
        const Real* const src_tracer, Real* const tgt_tracer, const Int ntracers,
        const Real* const src_density, Real* const tgt_density,
        const bool rho_isl) {
  assert(src_density && tgt_density);
  Timer::start(Timer::ts_isl_interp);
  const Int cnn = nslices(c2d.ptr) - 1;
  assert(cnn == nslices(advected_p));
  const Int dnn = nslices(m.dglln2cglln);
  const Int ne = std::sqrt(ncell_ / 6);
  const Int np = m.np, np2 = square(np);
  const bool footprint = ntracers > 0 && isl_impl_->footprint != nullptr;
  const Real* wt;
  auto& colscl = isl_impl_->wrk;
  if (rho_isl)
    calc_jacobian_departure(m, md_ != nullptr, advected_p, Jdep_.data());
  if (footprint) isl_impl_->footprint->start();
  ompparfor for (Int tni = 0; tni < cnn; ++tni) {
    // GLL values in the source cell corresponding to this departure point.
    Int ci;
    Real va[GLL::np_max], vb[GLL::np_max];
    Real a, b;
    calc_departure_data(m, *rd_, ne, tni, advected_p, ci, true, va, vb, a, b);
    cn_src_cell_[tni] = ci;
    if (rho_isl) {
      const Real* const src_cell = src_density + ci*np2;
      Real iv = 0;
      for (Int j = 0, k = 0; j < np; ++j)
        for (Int i = 0; i < np; ++i, ++k)
          iv += va[i]*vb[j]*src_cell[k];
      for (Int j = c2d.ptr[tni]; j < c2d.ptr[tni+1]; ++j) {
        const Int di = c2d.r[j];
        tgt_density[di] = (Jdep_(di)/Je_(di))*iv;
      }
    }
    for (Int tri = 0; tri < ntracers; ++tri) {
      // Interp the source cell to this point.
      const Real* const src_cell = src_tracer + tri*dnn + ci*np2;
      Real* const tgt = tgt_tracer + tri*dnn;
      Real iv = 0;
      for (Int j = 0, k = 0; j < np; ++j)
        for (Int i = 0; i < np; ++i, ++k)
          iv += va[i]*vb[j]*src_cell[k];
      // Scatter to DGLL target values corresponding to this CGLL node.
      for (Int j = c2d.ptr[tni]; j < c2d.ptr[tni+1]; ++j)
        tgt[c2d.r[j]] = iv;
      if (footprint && tri == 0) {
        for (Int j = c2d.ptr[tni]; j < c2d.ptr[tni+1]; ++j)
          isl_impl_->footprint->record(c2d.r[j], ci);
      }
    }
  }
  if (footprint) isl_impl_->footprint->finish();
  Timer::stop(Timer::ts_isl_interp);
}

void Remapper::
isl_cdr_rho (const Mesh& m, const RemapData& rd_src, const RemapData& rd_tgt,
             const MonoData& md, const spf::MassRedistributor::Ptr& mrd,
             Real* const src_density, const Int np2_src,
             Real* const tgt_density, const Int np2_tgt,
             const Int ntracers) {
  const Int len_src = np2_src*ncell_, len_tgt = np2_tgt*ncell_;
  Real* src, * tgt;
  dgll_set_src_tgt_in(src_density, tgt_density, src, tgt);
  const auto& F_src = rd_src.dgbfi_mass();
  const auto& F_tgt = rd_tgt.dgbfi_mass();
  Real mass_src = 0, mass_tgt = 0;
  accumulate_threaded_bfb<1>(
    [&] (const Int i, Real* a) { *a += F_src[i]*src[i]; },
    len_src, &mass_src);
  accumulate_threaded_bfb<1>(
    [&] (const Int i, Real* a) { *a += F_tgt[i]*tgt[i]; },
    len_tgt, &mass_tgt);
  if (Filter::is_nodal(md.global_type())) {
    glbl_only_pve(isl_impl_->gocd, len_tgt, F_tgt.data(),
                  mass_src - mass_tgt, tgt);
    dgll_set_tgt_out(tgt, tgt_density);
    if (record_total_mass_redistribution_)
      total_mass_redistribution_[ntracers] = total_mass_discrepancy_[ntracers] =
        std::abs(mass_src - mass_tgt);
    return;
  }
# pragma omp parallel for
  for (Int ti = 0; ti < ncell_; ++ti)
    mrd->record(ti, 0, 2, tgt, nullptr);
  // Global mass redistribution based on cell mean boundedness.
  mrd->redistribute_mass(mass_src - mass_tgt);
  if (record_total_mass_redistribution_) {
    Real tmr;
    accumulate_threaded_bfb<1>(
      [&] (const Int i, Real* a) { *a += std::abs(mrd->get_delta_mass(i)); },
      ncell_, &tmr);
    total_mass_redistribution_[ntracers] = tmr;
    total_mass_discrepancy_[ntracers] = std::abs(mass_src - mass_tgt);
  }
  // Shape preservation filter with cell mass deltas.
# pragma omp parallel for
  for (Int ti = 0; ti < ncell_; ++ti) {
    const Real rho_extra = mrd->get_delta_mass(ti);
    //Int info =
    md.limit_density(ti, tgt + ti*np2_tgt, rho_extra);
    //assert(info >= 0);
  }
  dgll_set_tgt_out(tgt, tgt_density);
}

// Constrained density reconstruction for classical semi-Lagrangian.
void Remapper::
isl_cdr (const Mesh& m, const RemapData& rd_src, const RemapData& rd_tgt,
         const MonoData& md, const spf::MassRedistributor::Ptr& mrd,
         Real* const src_density, Real* const src_tracer, const Int np_src,
         Real* const tgt_density, Real* const tgt_tracer, const Int ntracers,
         const bool positive_only, const bool rho_isl, const Filter::Enum cdr_method) {
  const auto& td = *isl_impl_->tend;
  const Int np_tgt = m.np, np2_tgt = square(np_tgt);
  const Int np2_src = square(np_src);
  const Int len_src = np2_src*ncell_, len_tgt = np2_tgt*ncell_;
  const auto& F_src = rd_src.dgbfi_mass();
  const auto& F_tgt = rd_tgt.dgbfi_mass();
  assert(Filter::is_nodal(md.global_type()) || mrd->get_np2() == np2_tgt);
  if (record_total_mass_redistribution_) {
    total_mass_redistribution_.resize(ntracers+1);
    total_mass_discrepancy_.resize(ntracers+1);
  }
  if (rho_isl) {
    isl_cdr_rho(m, rd_src, rd_tgt, md, mrd,
                src_density, np2_src,
                tgt_density, np2_tgt, ntracers);
    assert(nslices(m.dglln2cglln) == d2cer_->get_dnn());
    dss_rho(*d2cer_, tgt_density, isl_impl_->wrk);
  }
  if ( ! src_tracer) return;
  for (Int tri = 0; tri < ntracers; ++tri) {
    Real* src, * tgt;
    dgll_set_src_tgt_in(src_tracer + tri*len_src,
                        tgt_tracer + tri*len_tgt, src, tgt);
    // q -> Q.
    Real Q_mass_src, Q_mass_tgt;
    accumulate_threaded_bfb<1>(
      [&] (const Int i, Real* a) { *a += F_src[i]*src[i]*src_density[i]; },
      len_src, &Q_mass_src);
    accumulate_threaded_bfb<1>(
      [&] (const Int i, Real* a) {
        tgt[i] *= tgt_density[i];
        *a += F_tgt[i]*tgt[i];
      },
      len_tgt, &Q_mass_tgt);
    if (positive_only) {
      assert((md.local_type() != Limiter::none));
#     pragma omp parallel for
      for (Int ti = 0; ti < ncell_; ++ti)
        mrd->record(ti, 0, 2, tgt_density, tgt);
    } else {
      // Compute bounds for each source cell or subcell.
      SIQK_THROW_IF(subcell_bounds_, "not impl'ed");
#     pragma omp parallel for
      for (Int ci = 0; ci < ncell_; ++ci) {
        Real q_min, q_max;
        const Int os = ci*np2_src;
        q_min = q_max = src[os];
        for (Int i = 1; i < np2_src; ++i) {
          q_min = std::min(q_min, src[os + i]);
          q_max = std::max(q_max, src[os + i]);
        }
        if ( ! td.cdr_on_src_tends &&
            (td.dss_source || (tri >= td.idx_beg && tri < td.idx_end))) {
          // Get the tightest bounds we can.
          q_min = std::max(q_min, td.q_min[ncell_*tri + ci]);
          q_max = std::min(q_max, td.q_max[ncell_*tri + ci]);
        }
        if (fit_extremum_) {
          assert(np_tgt == np_); // for now
          const Int os = ci*np2_src;
          bool use;
          Real min, max;
          fit_extremum_->calc(src + os, min, max, use);
          if (use) {
            q_min = std::min(q_min, min);
            q_max = std::max(q_max, max);
          }
        }
        q_data_(ci,0) = q_min;
        q_data_(ci,1) = q_max;
      }
      // Compute bounds for each target node.
      if (mrd) {
#       pragma omp parallel for
        for (Int ti = 0; ti < ncell_; ++ti) {
          Real q_min[GLL::np_max*GLL::np_max], q_max[GLL::np_max*GLL::np_max];
          for (Int i = 0; i < np2_tgt; ++i) {
            const Int si = cn_src_cell_[m.dglln2cglln[ti*np2_tgt + i]];
            q_min[i] = q_data_(si,0);
            q_max[i] = q_data_(si,1);
          }
          mrd->record(ti, q_min, q_max, tgt_density, tgt);
        }
      }
    }
    const bool point = Filter::is_nodal(md.global_type());
    if ( ! point) {
      mrd->redistribute_mass(Q_mass_src - Q_mass_tgt);
      if (record_total_mass_redistribution_) {
        Real tmr;
        accumulate_threaded_bfb<1>(
          [&] (const Int i, Real* a) { *a += std::abs(mrd->get_delta_mass(i)); },
          ncell_, &tmr);
        total_mass_redistribution_[tri] = tmr;
        total_mass_discrepancy_[tri] = std::abs(Q_mass_src - Q_mass_tgt);
      }
    }
    // Shape preservation filter with cell mass deltas.
#   pragma omp parallel for
    for (Int ti = 0; ti < ncell_; ++ti) {
      Real* const tgti = tgt + ti*np2_tgt;
      const Real* const tgtdi = tgt_density + ti*np2_tgt;
      const Real Qm_extra = point ? 0 : mrd->get_delta_mass(ti);
      Real q_min[GLL::np_max*GLL::np_max], q_max[GLL::np_max*GLL::np_max];
      if (positive_only)
        md.limit_density(ti, tgti, Qm_extra);
      else {
        for (Int i = 0; i < np2_tgt; ++i) {
          const Int si = cn_src_cell_[m.dglln2cglln[ti*np2_tgt + i]];
          q_min[i] = q_data_(si,0);
          q_max[i] = q_data_(si,1);
          if (point) {
            const Real del = TendencyData::lcl_cdr_relax_fac*(q_max[i] - q_min[i]);
            q_min[i] -= del;
            q_max[i] += del;
          }
        }
        md.limit_tracer(ti, q_min, q_max, tgtdi, tgti, point /* expand */,
                        Qm_extra, true /* nocheck */);
      }
    }
    if (point) {
      Real extra_mass = Q_mass_src - Q_mass_tgt;
      glbl_only_lcldyn(isl_impl_->gocd, len_tgt, ncell_, q_data_.data(),
                       cn_src_cell_.data(), m.dglln2cglln.data(),
                       F_tgt.data(), tgt_density, extra_mass, tgt);
      if (record_total_mass_redistribution_) {
        total_mass_redistribution_[tri] = std::abs(extra_mass);
        total_mass_discrepancy_[tri] = std::abs(Q_mass_src - Q_mass_tgt);
      }
    }
#   pragma omp parallel for
    for (Int ti = 0; ti < ncell_; ++ti) {
      Real* const tgti = tgt + ti*np2_tgt;
      const Real* const tgtdi = tgt_density + ti*np2_tgt;
      Real q_min[GLL::np_max*GLL::np_max], q_max[GLL::np_max*GLL::np_max];
      if ( ! positive_only) {
        for (Int i = 0; i < np2_tgt; ++i) {
          const Int si = cn_src_cell_[m.dglln2cglln[ti*np2_tgt + i]];
          q_min[i] = q_data_(si,0);
          q_max[i] = q_data_(si,1);
        }
      }
      // Q -> q.
      for (Int i = 0; i < np2_tgt; ++i) {
        if (tgtdi[i] == 0) tgti[i] = positive_only ? 0 : q_min[i];
        else {
          tgti[i] /= tgtdi[i];
          if ( ! positive_only) {
            // Clean up numerical error.
            tgti[i] = std::max(q_min[i], std::min(q_max[i], tgti[i]));
          }
        }
      }
    }
    dgll_set_tgt_out(tgt, tgt_tracer + tri*len_tgt);    
  }
}

void Remapper
::init_physgrid (const Real* const rho, const Real* const tracer, const Int ntracers,
                 pg::PhysgridData& pg) {
  const auto cdr_method = md_ ? md_->global_type() : Filter::none;
  const bool run_cdr = cdr_method != Filter::none;
  if (pr_->experiment <= 1) {
    isl_impl_->transfer_to_physgrid(pg, run_cdr, ntracers, rho, tracer);
  } else {
    isl_impl_->init_tracer_arrays(ntracers);
    Real* const rho_impl = isl_impl_->rho[isl_impl_->idxs[0]].data();
    Real* const tracer_impl = isl_impl_->tracer[isl_impl_->idxs[0]].data();
    isl_impl_->interp_rho(rho, rho_impl, false);
    isl_impl_->transfer_q_to_t_mesh(run_cdr, ntracers,
                                    rho, tracer,
                                    rho_impl, tracer_impl);
    isl_impl_->transfer_to_physgrid(pg, run_cdr, ntracers,
                                    rho_impl, tracer_impl);
  }  
}

void Remapper
::isl (const AVec3s& advected_p, Real* const src_tracer,
       Real* const tgt_tracer, const Int ntracers,
       Real* const src_density, Real* const tgt_density,
       const bool positive_only, const bool rho_isl,
       const Real* srcterm_tendency /* mixing ratio */,
       const Int srcterm_idx_beg, const Int srcterm_idx_end,
       const pg::PhysgridData::Ptr& pg) {
  Timer::start(Timer::ts_remap_project);
  const bool first = isl_impl_->start_time_step(
    ntracers, srcterm_tendency, srcterm_idx_beg, srcterm_idx_end, pg);
  const auto cdr_method = md_ ? md_->global_type() : Filter::none;
  const bool run_cdr = cdr_method != Filter::none;
  assert( ! (run_cdr && ! md_));
  const bool cdr_glbl_only = run_cdr && Filter::is_nodal(md_->global_type());
  // We'd need continuous p-refined rho if we were to run params on the
  // p-refined GLL grid, but in practice we're not interested in that; we want
  // to use physgrid. Set this to false so that I don't accidentally force rho
  // to be continuous. I prefer the simplicity of not doing that.
  const bool rho_continuous = false; //! pg;
  if (run_cdr && ! cdr_glbl_only && ! mrd_) {
    const auto& F = rd_->dgbfi_mass();
    const Int ne = std::sqrt(ncell_ / 6);
    mrd_ = std::make_shared<spf::MassRedistributor>(
      ne, F, Filter::to_mrd(cdr_method));
  }
  if (pg && pr_->experiment <= 1) {
    // Map physgrid tendencies to v-grid tracers.
    isl_impl_->add_tendencies_from_physgrid(
      src_density, src_tracer, run_cdr, ntracers, *isl_impl_->tend,
      isl_impl_->tend->dq, *pg, d2cer_);
    isl_impl_->dss_q(src_density, src_tracer, isl_impl_->tend->idx_beg,
                     isl_impl_->tend->idx_end);
  }
  if (pr_->experiment == 0) {
    // Conventional case, prior to any p-refinement experiments. Also does
    // original basic p-refinement.
    interp(*m_, c2d_v_, advected_p,
           src_tracer, tgt_tracer, ntracers,
           src_density, tgt_density, rho_isl);
    if (run_cdr)
      isl_cdr(*m_, *rd_, *rd_, *md_, mrd_,
              src_density, src_tracer, np_,
              tgt_density, tgt_tracer, ntracers,
              positive_only, rho_isl, cdr_method);
    else
      dss_rho(*d2cer_, tgt_density, isl_impl_->wrk);
    if (pg)
      isl_impl_->transfer_to_physgrid(*pg, run_cdr, ntracers, tgt_density, tgt_tracer);
  } else if (pr_->experiment == 1) {
    // slmmir sees the np>4 grid as primary, e.g. for error calculations and toy
    // chemistry. But, here, integrate rho on the v-grid and interpolate it to
    // the fine grid. Interpolation is done such that it's element local but
    // gives the same element boundary values.
    Real* const src_rho_impl = isl_impl_->rho[isl_impl_->idxs[0]].data();
    Real* const tgt_rho_impl = isl_impl_->rho[isl_impl_->idxs[1]].data();
    if (first) {
      Int n = ncell_ * square(pr_->m_f->np);
      for (Int i = 0; i < n; ++i)
        SIQK_THROW_IF(src_density[i] != 1,
                      "rho should be uniformly 1 in the first time step.");
      n = ncell_ * square(pr_->m_v->np);
      for (Int i = 0; i < n; ++i) src_rho_impl[i] = 1;
      if (run_cdr && ! isl_impl_->mrd_c && ! cdr_glbl_only) {
        const auto& F = pr_->rd_v->dgbfi_mass();
        const Int ne = std::sqrt(ncell_ / 6);
        isl_impl_->mrd_c = std::make_shared<spf::MassRedistributor>(
          ne, F, Filter::to_mrd(cdr_method));
      }
      d2cer_ = std::make_shared<D2Cer>(pr_->m_v->dglln2cglln, pr_->rd_v->dgbfi_mass());
    }
    if (rho_isl) {
      // Transport rho on v mesh to mimic dycore.
      interp(*pr_->m_v, c2d_v_, advected_p,
             nullptr, nullptr, 0,
             src_rho_impl, tgt_rho_impl, rho_isl);
      if (run_cdr) {
        // Property-preserve rho on v mesh to mimic dycore.
        isl_cdr(*pr_->m_v, *pr_->rd_v, *pr_->rd_v, *pr_->md_v, isl_impl_->mrd_c,
                src_rho_impl, nullptr, pr_->m_v->np,
                tgt_rho_impl, nullptr, ntracers,
                positive_only, rho_isl, cdr_method);
      } else {
        dss_rho(*d2cer_, tgt_rho_impl, isl_impl_->wrk);
      }
      // Interpolate current rho to f mesh.
      isl_impl_->interp_rho(tgt_rho_impl, tgt_density, rho_continuous);
    }  
    // Advect tracers on fine mesh.
    const auto ap_fine = isl_impl_->interp_p(advected_p);
    interp(*m_, c2d_f_, ap_fine,
           src_tracer, tgt_tracer, ntracers,
           // These are unused unless we're using the local mass conservation
           // option.
           src_rho_impl, tgt_rho_impl, false);
    if (run_cdr)
      isl_cdr(*m_, *rd_, *rd_, *md_, mrd_,
              src_density, src_tracer, np_,
              tgt_density, tgt_tracer, ntracers,
              positive_only, false, cdr_method);
    if (pg)
      isl_impl_->transfer_to_physgrid(*pg, run_cdr, ntracers, tgt_density, tgt_tracer);
  } else if (pr_->experiment == 5) {
    // The v-grid is primary. Fine-grid data are internal to this impl.
    Real* const src_rho_impl = isl_impl_->rho[isl_impl_->idxs[0]].data();
    Real* const tgt_rho_impl = isl_impl_->rho[isl_impl_->idxs[1]].data();
    Real* const src_tracer_impl = isl_impl_->tracer[isl_impl_->idxs[0]].data();
    Real* const tgt_tracer_impl = isl_impl_->tracer[isl_impl_->idxs[1]].data();
    if (rho_isl) {
      // Transport rho on v mesh to mimic dycore.
      interp(*pr_->m_v, c2d_v_, advected_p,
             nullptr, nullptr, 0,
             src_density, tgt_density, rho_isl);
      if (run_cdr) {
        // Property-preserve rho on v mesh to mimic dycore.
        isl_cdr(*pr_->m_v, *pr_->rd_v, *pr_->rd_v, *pr_->md_v, mrd_,
                src_density, nullptr, pr_->m_v->np,
                tgt_density, nullptr, ntracers,
                positive_only, rho_isl, cdr_method);
      } else {
        dss_rho(*d2cer_, tgt_density, isl_impl_->wrk);
      }
    }
    if (first) {
      isl_impl_->interp_rho(src_density, src_rho_impl, rho_continuous);
      isl_impl_->transfer_q_to_t_mesh(run_cdr, ntracers,
                                      src_density, src_tracer,
                                      src_rho_impl, src_tracer_impl);
      if (run_cdr && ! isl_impl_->mrd_f && ! cdr_glbl_only) {
        const auto& F = pr_->rd_f->dgbfi_mass();
        const Int ne = std::sqrt(ncell_ / 6);
        isl_impl_->mrd_f = std::make_shared<spf::MassRedistributor>(
          ne, F, Filter::to_mrd(cdr_method));
      }
    }
    // Map tendency from physics (pg or v) to t grid.
    isl_impl_->add_tendencies(src_density, src_tracer,
                              src_rho_impl, src_tracer_impl,
                              run_cdr, ntracers, *isl_impl_->tend,
                              isl_impl_->tend->dq, pg);
    // Either or both of pg->tgrid source terms or caas-node with relaxed lcl
    // caas prefilter causes q to be discontinuous at the end of the advection
    // step. Restore continuity.
    isl_impl_->dss_q(src_rho_impl, src_tracer_impl, 0, ntracers);
    // Advect tracers on fine mesh.
    const auto ap_fine = isl_impl_->interp_p(advected_p);
    interp(*pr_->m_f, c2d_f_, ap_fine,
           src_tracer_impl, tgt_tracer_impl, ntracers,
           // These are unused unless we're using the local mass conservation
           // option.
           src_rho_impl, tgt_rho_impl, false);
    // Interpolate current rho to f mesh.
    isl_impl_->interp_rho(tgt_density, tgt_rho_impl, rho_continuous);
    // Global property preservation on t mesh, local on both.
    if (run_cdr) {
      // Property preserve q on f mesh.
      isl_cdr(*pr_->m_f, *pr_->rd_f, *pr_->rd_f, *pr_->md_f, isl_impl_->mrd_f,
              src_rho_impl, src_tracer_impl, pr_->m_f->np,
              tgt_rho_impl, tgt_tracer_impl, ntracers,
              positive_only, false /* don't apply cdr to rho */, cdr_method);
    }
    // In the following two cases, continuity does not hold nor need hold on
    // any of the grids.
    //   Map q on tgrid to q on vgrid.
    isl_impl_->transfer_q_to_v_mesh(run_cdr, ntracers,
                                    tgt_rho_impl, tgt_tracer_impl,
                                    tgt_density, tgt_tracer);
    //   Map q on tgrid to q on pgrid.
    if (pg)
      isl_impl_->transfer_to_physgrid(*pg, run_cdr, ntracers,
                                      tgt_rho_impl, tgt_tracer_impl);
  } else {
    SIQK_THROW_IF(true, "not supported; wasn't useful");
  }
  isl_impl_->finish_time_step();
  Timer::stop(Timer::ts_remap_project);
}

void Remapper::
get_prefine_internal (const Real*& q, const Real*& dgbfi_mass, Int& len,
                      const bool first) const {
  q = isl_impl_->tracer[isl_impl_->idxs[first ? 1 : 0]].data();
  dgbfi_mass = pr_->rd_f->dgbfi_mass().data();
  len = pr_->rd_f->dgbfi_mass().size();
}
