#include "mpas_advect_ir_test.hpp"

#include <set>
#include <iostream>
#include <sstream>

namespace ir {

const InitialCondition::Shape TestProblemOptions::all_ics[max_n_ics] = {
  InitialCondition::Shape::XYZTrig, InitialCondition::Shape::GaussianHills,
  InitialCondition::Shape::CosineBells, InitialCondition::Shape::SlottedCylinders
};

std::string to_string (const TestProblemOptions& tpo, const Params& p,
                       const Comm& he) {
  std::stringstream ss;
  ss << "(" << p.vf->name() << " "
     << (tpo.ics_type == IcsType::all_2d ? "all" :
         tpo.ics_type == IcsType::dcmip12t11 ? "dcmip12t11" :
         InitialCondition::to_string(tpo.ic)) << " "
     << Trajectory::Method::to_string(p.advect_method) << " nstep " << p.nstep
     << " nsubstep " << p.nsubstep << " trajhalo " << p.traj_nhalo
     << " commhalo " << he.szhalo() << " order " << p.order << " tq_order "
     << p.tq_order << " limit " << p.limit << " cedr " << p.cedr << " perturb "
     << p.perturb_dep_pts << " ll_cell_int " << p.lag_lvls_cell_int
     << " bfb " << p.bfb << ")";
  return ss.str();
}

void optstr2options (std::string optstr, TestProblemOptions& tpo, Params& o,
                     const bool verbose=true) {
  tpo.ic = InitialCondition::Shape::XYZTrig;
  tpo.ics_type = IcsType::one_2d;
  if (optstr.size() == 0) optstr = "r";
  o.vf = vf_from_short(optstr[0]);
  do {
    if (optstr.size() == 1) break;
    if (optstr[1] == 'a') tpo.ics_type = IcsType::all_2d;
    else if (optstr[1] == 'A') tpo.ics_type = IcsType::dcmip12t11;
    else tpo.ic = ic_from_short(optstr[1], verbose);
    if (optstr.size() == 2) break;
    o.order = optstr[2] - '0';
    if (optstr.size() == 3) break;
    o.advect_method = Trajectory::Method::from_short(optstr[3]);
    if (optstr.size() == 4) break;
    o.limit = optstr[4] - '0';
    if (optstr.size() == 5) break;
    o.cedr = optstr[5] - '0';
    if (optstr.size() == 6) break;
    o.perturb_dep_pts = optstr[6] - '0';
  } while (0);
  if (o.order < 1 || o.order > 3) {
    if (verbose)
      printf("slot 2: order must be in 1,2,3; is %d; setting to 2\n", o.order);
    o.order = 2;
  }
  if (o.limit && o.order != 2) {
    if (verbose)
      printf("slot 5: limit works only with order 2; setting to off\n");
    o.limit = 0;
  }
}

static void fill_ics (const MpasMesh& mm, const TestProblemOptions& tpo,
                      const Params& p, Ar2& rho, Ar3& q) {
  assert(mm.l.get_interface_z(0));
  const int nlvl = mm.l.nmax;
  ompparfor for (Int ci = 0; ci < mm.t.ncell; ++ci) {
    const Real* zint = mm.l.get_interface_z(ci);
    for (Int li = 0; li < nlvl; ++li) {      
      const Real zeta0 = zint[li], zeta1 = zint[li+1];
      rho(li,ci) = p.vf->get_rho(zeta0, zeta1);
      const Real profile = rho(li,ci)*p.vf->get_ic_shape_fn(zeta0, zeta1);
      switch (tpo.ics_type) {
      case IcsType::one_2d:
        q(0,li,ci) = profile*get_ic(mm, tpo.ic, ci);
        break;
      case IcsType::all_2d:
        for (int qi = 0; qi < p.ntracer; ++qi)
          q(qi,li,ci) = profile*get_ic(mm, TestProblemOptions::all_ics[qi], ci);
        break;
      case IcsType::dcmip12t11:
        for (int qi = 0; qi < p.ntracer; ++qi)
          q(qi,li,ci) = profile*get_ic_dcmip12t11(mm, qi, ci, zeta0, zeta1);
        break;
      }
    }
  }
}

void set (const MpasMesh& mm, Comm& he, TestProblemOptions& o, Params& p) {
  optstr2options(o.flow_ic_str, o, p, he.amroot());
  p.ntracer = IcsType::get_ntracer(o.ics_type);

  p.nstep = o.nstep;
  if (p.nstep < 1) p.nstep = calc_default_nstep(mm, he, *p.vf);
  p.nsubstep = o.nsubstep;
  p.traj_nhalo = o.traj_nhalo;
  calc_substep_data(mm, he, *p.vf, p.nstep, p.traj_nhalo, p.nsubstep);

  if (he.amroot()) printf("optstr %s\n", o.flow_ic_str.c_str());
}

static void write (const OutputMgr& om, const Real time, const Acr2& rho,
                   const Acr3& Q) {
  om.ncw->advance_time_to(time);
  for (const auto& d : om.q_descs) om.ncw->write_field(d, &rho, Q);
  for (const auto& d : om.rho_descs) om.ncw->write_field(d, rho);
}

#if 0
static void summarize_rho (Comm& he, const Real t, const int ncell, const Acr2& rho) {
  Real lmin = 1e10, lmax = -1e10;
  for (int k = 0; k < rho.m(); ++k)
    for (int i = 0; i < ncell; ++i) {
      lmin = std::min(lmin, rho(k,i));
      lmax = std::max(lmax, rho(k,i));
    }
  Real gmin, gmax;
  he.all_reduce(Comm::MIN, 1, &lmin, &gmin);
  he.all_reduce(Comm::MAX, 1, &lmax, &gmax);
  if (he.amroot())
    printf("day %5.2f min %6.3f max %6.3f\n", t/(24*3600), gmin, gmax);
}
#endif

void run_test_problem (const RemapFormOp::Ptr& op, TestProblemOptions& o,
                       const OutputMgr::ConstPtr om, Real* ic, Real* mid, Real* fin) {
  Timer::init();
  Timer::start(Timer::total);

  const auto mm = op->get_mpas_mesh();
  const auto he = op->get_comm();
  const auto& p = op->get_params();
  const int ncell = mm->t.ncell;

  if (he->amroot()) std::cout << to_string(o, p, *he) << "\n";

  const Int nmlvl = mm->l.nmax;
  Ar2 rhos[2];
  for (int i = 0; i < 2; ++i) rhos[i].reset(nmlvl, ncell);
  Ar3 qs[2];
  for (int i = 0; i < 2; ++i) qs[i].reset(p.ntracer, nmlvl, ncell);

  fill_ics(*mm, o, p, rhos[0], qs[0]);
  if (om) write(*om, 0, rhos[0], qs[0]);

  const bool measure_q = ! (p.perturb_dep_pts && ! p.cedr);

  if (ic) {
    auto& q = qs[0];
    auto& rho = rhos[0];
    for (int i = 0; i < ncell; ++i) ic[i] = q(0,0,i)/(measure_q ? rho(0,i) : 1);
  }

  const Real te = p.vf->cycle_time(), dt = te/p.nstep;
  int i0 = 0, i1 = 1;
  for (int ti = 0; ti < p.nstep; ++ti) {
    const Real ts = (Real(ti)/p.nstep)*te;
    //summarize_rho(*he, ts, mm->ncown, rhos[i0]);
    op->integrate(ts, dt, rhos[i0], rhos[i1], qs[i0], qs[i1]);
    auto& rho = rhos[i1];
    auto& q = qs[i1];
    if (mid && ti == p.nstep/2-1)
      for (int i = 0; i < ncell; ++i)
        mid[i] = q(0,0,i)/(measure_q ? rho(0,i) : 1);
    if (o.movie)
      for (int i = 0; i < ncell; ++i)
        ic[ncell*(ti+1) + i] = q(0,0,i)/(measure_q ? rho(0,i) : 1);
    if (om && ((ti+1) % om->time_step_write_stride == 0 || ti == p.nstep-1))
      write(*om, (Real(ti+1)/p.nstep)*te, rho, q);
    std::swap(i0, i1);
  }

  if (fin) {
    auto& q = qs[i0];
    auto& rho = rhos[i0];
    for (int i = 0; i < ncell; ++i) fin[i] = q(0,0,i)/(measure_q ? rho(0,i) : 1);
  }

  fill_ics(*mm, o, p, rhos[i1], qs[i1]);

  for (const auto cmp_q : {true, false}) {
    if ( ! measure_q && cmp_q) continue;
    CmpOutput e;
    std::string qtag(cmp_q ? "q" : "Q");
    const auto f = [&] (const int li) {
      for (Int qi = 0; qi < p.ntracer; ++qi) {
        cmp(*mm, *he, qi, li, nmlvl, mm->ncown,
            cmp_q ? &rhos[i1] : nullptr, qs[i1],
            cmp_q ? &rhos[i0] : nullptr, qs[i0],
            e, mm->g.my_areaCell.data());
        if (he->amroot()) {
          std::string mass_msg, lim_msg;
          if (p.vf->is_3d() && ! p.cedr) mass_msg = " (3D and no CEDR)";
          if ( ! (p.limit || p.cedr)) lim_msg = " (no limiter)";
          printf("re %d %s %d l2 %1.16e l1 %1.16e linf %1.16e\n"
                 "re mass %1.16e%s\n"
                 "lim %1.16e %1.16e%s\n",
                 li, qtag.c_str(), qi, e.l2, e.l1, e.linf,
                 e.mass, mass_msg.c_str(),
                 e.lo, e.hi, lim_msg.c_str());
        }
      }      
    };
    if (p.vf->is_3d())
      f(-1);
    else {
      for (Int li = 0; li < nmlvl; ++li)
        f(li);
    }
  }

  Timer::stop(Timer::total);
  Timer::print();
}

bool unittest (const MpasMesh::ConstPtr& mm, const siqp::sh::Mesh::ConstPtr& m,
               const Comm::Ptr& he, Real* data) {
  Params p;
  p.nstep = 8;
  p.nsubstep = 6;
  p.traj_nhalo = 6;
  p.ntracer = 1;
  TestProblemOptions o;
  optstr2options("nt22", o, p, he->amroot());
  const auto op = std::make_shared<RemapFormOp>(mm, m, he, p);
  return op->unittest(data);
}

void advect_mesh (const MpasMesh::ConstPtr& mm, const siqp::sh::Mesh::ConstPtr& m,
                  TestProblemOptions& o,
                  const int stepno, std::vector<Real>& xyzs_orig,
                  std::vector<Real>& xyzs_advected) {
  Params p;
  p.nstep = o.nstep;
  optstr2options(o.flow_ic_str, o, p);
  assert(stepno <= p.nstep);

  const auto he = Comm::create();

  const Real te = p.vf->cycle_time(), dt = te/p.nstep;
  p.nsubstep = -1;
  p.traj_nhalo = -1;
  calc_substep_data(*mm, *he, *p.vf, p.nstep, p.traj_nhalo, p.nsubstep);

  const auto op = std::make_shared<RemapFormOp>(mm, m, he, p);

  op->advect_mesh((Real(stepno-1)/p.nstep)*te, dt);
  op->fill_patches(xyzs_orig, xyzs_advected);
}

} // namespace ir
