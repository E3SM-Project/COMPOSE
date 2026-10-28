#include "mpas_lagrangian_levels.hpp"
#include "mpas_gallery.hpp"
#include "mpas_geometry.hpp"
#include "mpas_vr.hpp"
#include "mpas_accum.hpp"

// Index names have the form
//     {(c)ell,(l)evel,(v)ertex}{(l)agrangian,(e)ulieran}{index letter}

struct LagrangianLevels::Impl {
  static constexpr int order = 2;
  // Inexact integration is fine.
  static constexpr int tq_order = 2;

  MpasMesh::ConstPtr mm;
  siqp::sh::Mesh::ConstPtr m;
  Comm::Ptr he;
  Trajectory::ConstPtr traj;
  OctreeMgr::ConstPtr otm;
  Params params;
  CellFromVtxsShapeFunction cvsf;
  Real local_accum_maxmag;

  void init(const MpasMesh::ConstPtr& mm, const siqp::sh::Mesh::ConstPtr& m,
            const Comm::Ptr& he, const Trajectory::ConstPtr& traj,
            const OctreeMgr::ConstPtr& otm, const Params& p);

  void transfer_lv2ec(const Acr3& plag, const Acr3& f_vtx, const Ar1& wrk,
                      const Ar2& f_cc);
  void remap_vertically(const Acr2& zeta_eul, const Ar2& zeta_lag, const Ar2& rho,
                        const Ar1& rwrk, const int nq, const Ar3& Q);

  // For testing.
  void transfer_ec2ev(const CellFromCellsModalShapeFunction& ccsf,
                      const Acr2& f_cc, const Ar3& f_vtx);

private:
  void search(const Acr3& plag, const int cli, const int lvl, OTSearchFunctor& f);
  template <typename T> void accumulate_into_ec(
    const Acr3& plag, const Acr3& f_vtx, const Array2D<T>& f_cc);
  void interp_to_ec(const Acr3& plag, const Acr3& f_vtx, const Ar2& f_cc);
};

static void get_max_3d_params (const MpasMesh& mm, Comm& he,
                               Real& max_cell_area, Real& max_zeta) {
  const auto& geo = mm.g;
  const auto& lvl = mm.l;
  Real l[2] = {0}, g[2];
  for (int ci = 0; ci < mm.ncown; ++ci) {
    l[0] = std::max(l[0], geo.my_areaCell[ci]);
    l[1] = std::max(l[1], lvl.get_interface_z(ci)[lvl.get_cmax(ci)+1]);
  }
  he.all_reduce(Comm::MAX, 2, l, g);
  max_cell_area = g[0];
  max_zeta = g[1];
}

void LagrangianLevels::Impl
::init (const MpasMesh::ConstPtr& mm_, const siqp::sh::Mesh::ConstPtr& m_,
        const Comm::Ptr& he_, const Trajectory::ConstPtr& traj_,
        const OctreeMgr::ConstPtr& otm_, const Params& params_) {
  mm = mm_;
  m = m_;
  he = he_;
  traj = traj_;
  otm = otm_;
  params = params_;
  cvsf.init(mm, he, order, tq_order, true);

  if (params.cell_integrated && params.bfb) {
    Real max_cell_area, max_zeta;
    get_max_3d_params(*mm, *he, max_cell_area, max_zeta);
    //todo Improve numerics. This is the coarsest, most straightforward bound.
    local_accum_maxmag = max_cell_area*max_zeta;
  }
}

void LagrangianLevels::Impl
::transfer_ec2ev (const CellFromCellsModalShapeFunction& ccsf, const Acr2& f_cc,
                  const Ar3& f_vtx) {
  const auto ncown = mm->ncown, nlvl = traj->num_levels();
  const Acr3 f_cc_3(1, f_cc.m(), f_cc.n(), f_cc.data());
  ompparfor for (Int ci = 0; ci < ncown; ++ci)
    for (Int li = 0; li < nlvl; ++li) {
      Real shape[CellShapeFunction::max_ho_ncoef];
      ccsf.calc_shape(ci, f_cc_3, 0, li, shape);
      ccsf.eval_shape_at_vtxs(ci, f_cc_3(0,li,ci), shape, 0, li, f_vtx(li,ci));
    }
}

// Collect potential overlapping Eulerian cells.
void LagrangianLevels::Impl
::search (const Acr3& plag, const int cli, const int lvl, OTSearchFunctor& f) {
  const auto& t = mm->t;
  const auto& e = m->e;
  const int nvl = t.c2ne[cli];
  if (traj->has_in_cell_data()) {
    for (Int vli = 0; vli < nvl; ++vli) {
      const auto cei = traj->
        get_cell_lvl_containing_vtx(t.get_c2v(cli,vli), lvl).cell();
      f(cei);
      const auto nee = t.c2ne[cei];
      for (Int i = 0; i < nee; ++i)
        f(t.get_c2c(cei,i));
    }
  } else {
    otm->apply(plag, lvl, e(cli), e.m(), f);
  }
  require_msg(f.ok, "LagrangianLevels::search: Too many hits.");
  f.post();
}

template <typename T> inline void
accum_integral (const Array2D<T>& f_cc, const int cei, const int lvl,
                const Real integral) {
  auto& f_cc_lvl_cei = f_cc(lvl,cei);
  f_cc_lvl_cei += integral;
}

static inline void
accum_integral (const Ar2& f_cc, const int cei, const int lvl,
                const Real integral) {
  auto& f_cc_lvl_cei = f_cc(lvl,cei);
  ompatomic f_cc_lvl_cei += integral;
}

// Accumulate Lag own src cell contributions to owned+haloed Eul tgt cells.
template <typename AccumT>
void LagrangianLevels::Impl
::accumulate_into_ec (const Acr3& plag, const Acr3& f_vtx,
                      const Array2D<AccumT>& f_cc) {
  const auto nlvl = traj->num_levels(), ncown = mm->ncown;
  const auto& t = mm->t;
  const auto& l = mm->l;
  const auto& e = m->e;
  const Acr4 f_vtx_4(f_vtx.m(), 1, f_vtx.n(), f_vtx.o(), f_vtx.data());
  TriangleQuadrature tq;

  for (int lvl = 0; lvl < nlvl; ++lvl) {
    cvsf.reinit(plag, lvl);

    ompparfor for (Int cli = 0; cli < ncown; ++cli) {
      if (lvl < l.get_cmin(cli) || lvl > l.get_cmax(cli) + 1) continue;

      OTSearchFunctor f;
      search(plag, cli, lvl, f);

      Real shape[CellShapeFunction::max_ho_ncoef];
      cvsf.calc_shape(cli, 0, f_vtx_4, 0, lvl, shape);

      const int nvl = t.c2ne[cli];
      PolyClipper pc;
      pc.ni = nvl;
      for (int i = 0; i < nvl; ++i)
        geo::copy(pc.vi(i), plag(lvl,e(i,cli)));

      for (Int hi = 0; hi < f.nhit; ++hi) {
        const auto cei = f.hits[hi];
        if (cei == -1) continue;
        pc.clip(*m, cei);
        if (pc.no == 0) continue;
        const Real integral = cvsf.integrate_over_polygon(tq, pc.vo, pc.no,
                                                          shape, cli);
        accum_integral(f_cc, cei, lvl, integral);
      }
    }
  }
}

// Interpolate Lag own src cell zeta to owned+haloed Eul tgt cells.
void LagrangianLevels::Impl
::interp_to_ec (const Acr3& plag, const Acr3& f_vtx, const Ar2& f_cc) {
  const auto nlvl = traj->num_levels(), ncown = mm->ncown;
  const auto& t = mm->t;
  const auto& l = mm->l;
  const auto& ccs = mm->g.my_cellCentroid;
  const auto& e = m->e;
  const Acr4 f_vtx_4(f_vtx.m(), 1, f_vtx.n(), f_vtx.o(), f_vtx.data());

  for (int lvl = 0; lvl < nlvl; ++lvl) {
    cvsf.reinit(plag, lvl);

    ompparfor for (Int cli = 0; cli < ncown; ++cli) {
      if (lvl < l.get_cmin(cli) || lvl > l.get_cmax(cli) + 1) continue;

      OTSearchFunctor f;
      search(plag, cli, lvl, f);

      Real shape[CellShapeFunction::max_ho_ncoef];
      cvsf.calc_shape(cli, 0, f_vtx_4, 0, lvl, shape);

      // Lagrangian cell edge inward normals.
      const Int nvl = t.c2ne[cli];
      const auto elem = e(cli);
      Real nmls_buf[3*max_vtx];
      Ar2 nmls(3, nvl, nmls_buf);
      for (Int i = 0; i < nvl; ++i)
        geo::edge_normal(plag(lvl,elem[i]), plag(lvl,elem[(i+1)%nvl]), nmls(i));

      for (Int hi = 0; hi < f.nhit; ++hi) {
        const auto cei = f.hits[hi];
        if (cei == -1) continue;
        const auto cec = ccs(cei);
        //todo-pebfb Tie breaking.
        //todo-fp Guard against no f_cc value.
        bool outside = false;
        for (Int i = 0; i < nvl; ++i)
          if ( ! geo::inside(cec, plag(lvl,elem[i]), nmls(i))) {
            outside = true;
            break;
          }
        if (outside) continue;
        f_cc(lvl,cei) = cvsf.eval_shape_fn(cli, shape, cec);
      }
    }
  }
}

void LagrangianLevels::Impl
::transfer_lv2ec (const Acr3& plag, const Acr3& f_vtx, const Ar1& wrk,
                  const Ar2& f_cc) {
  const auto nlvl = traj->num_levels(), ncell = mm->t.ncell, ncown = mm->ncown;
  const auto& g = mm->g;
  if (params.cell_integrated) {
    if (params.bfb) {
      Array2D<OrderInvariantAccumTraits::LocalAtomic>
        f_cc_accum(f_cc.m(), f_cc.n());
      ompparfor for (int ci = 0; ci < ncell; ++ci)
        for (int li = 0; li < nlvl; ++li) {
          f_cc_accum(li,ci).set_maxmag(local_accum_maxmag);
          f_cc_accum(li,ci).set_irep(0);
        }
      accumulate_into_ec(plag, f_vtx, f_cc_accum);
      he->exchange_and_accumulate_cell(nlvl, f_cc_accum);
      ompparfor for (int ci = 0; ci < ncell; ++ci)
        for (int li = 0; li < nlvl; ++li)
          f_cc(li,ci) = f_cc_accum(li,ci)();
    } else {
      // Zero owned+haloed Eulerian tgt cell data.
      ompparfor for (int ci = 0; ci < ncell; ++ci)
        for (int li = 0; li < nlvl; ++li)
          f_cc(li,ci) = 0;
      accumulate_into_ec(plag, f_vtx, f_cc);
      // Comm non-0 values.
      he->exchange_and_accumulate_cell(nlvl, f_cc);
    }
    // Divide by area.
    ompparfor for (int ci = 0; ci < ncown; ++ci)
      for (int li = 0; li < nlvl; ++li)
        f_cc(li,ci) /= g.my_areaCell[ci];
  } else {
    ompparfor for (int ci = 0; ci < ncell; ++ci)
      for (int li = 0; li < nlvl; ++li)
        f_cc(li,ci) = 0;
    interp_to_ec(plag, f_vtx, f_cc);
    // Comm non-0 values.
    he->exchange_and_accumulate_cell(nlvl, f_cc);
  }
}

void LagrangianLevels::Impl
::remap_vertically (const Acr2& zeta_eul, const Ar2& zeta_lag, const Ar2& rho,
                    const Ar1& wrk, const int nq, const Ar3& Q) {
#define forli for (int li = 0; li < nmlvl; ++li)
  const auto nilvl = traj->num_levels(), nmlvl = nilvl-1, ncown = mm->ncown;
  assert(wrk.m() >= nmlvl*4*mpas_omp_get_max_threads());
  Ar3 rwrk(nmlvl, 4, mpas_omp_get_max_threads(), wrk.data());  
  bool ok = true;
  ompparfor for (int ci = 0; ci < ncown; ++ci) {
#ifndef NDEBUG
    if (ok) {
      ok = vr::assert_grids(nilvl, zeta_eul(ci), zeta_lag(ci), true);
      if ( ! ok) {
        prarr("zeul",zeta_eul(ci),nilvl);
        prarr("zlag",zeta_lag(ci),nilvl);
      }
    }
#endif
    const int tid = mpas_omp_get_thread_num();
    Real* const dps = rwrk(0,tid);
    Real* const dpt = rwrk(1,tid);
    Real* const fdp = rwrk(2,tid);
    Real* const rho_orig = rwrk(3,tid);
    // Remap rho. In practice rho is given to us.
    forli dps[li] = zeta_eul(li+1,ci) - zeta_eul(li,ci);
    forli dpt[li] = zeta_lag(li+1,ci) - zeta_lag(li,ci);
    forli rho_orig[li] = rho(li,ci);
    forli fdp[li] = rho_orig[li]*dps[li];
    vr::homme_vr(nmlvl, dps, dpt, fdp, false);
    forli rho(li,ci) = fdp[li]/dpt[li];
    // Now include rho in dp for Q.
    forli dpt[li] = fdp[li];
    forli dps[li] *= rho_orig[li];
    for (int qi = 0; qi < nq; ++qi) {
      forli fdp[li] = Q(qi,li,ci)*(zeta_eul(li+1,ci) - zeta_eul(li,ci));
      vr::homme_vr(nmlvl, dps, dpt, fdp, false);
#ifndef NDEBUG
      Real smass = 0, tmass = 0;
      forli smass += Q(qi,li,ci)*(zeta_eul(li+1,ci) - zeta_eul(li,ci));
      forli tmass += fdp[li];
      const Real re = reldif(smass, tmass);
      if (re > 1e-10) {
        pr(puf(ci) pu(qi) pu(smass) pu(tmass) pu(re));
        ok = false;
      }
#endif
      forli Q(qi,li,ci) = fdp[li]/(zeta_lag(li+1,ci) - zeta_lag(li,ci));
    }
  }
  assert(ok);
#undef forli
}

LagrangianLevels
::LagrangianLevels (const MpasMesh::ConstPtr& mm, const siqp::sh::Mesh::ConstPtr& m,
                    const Comm::Ptr& he, const Trajectory::ConstPtr& traj,
                    const OctreeMgr::ConstPtr& otm, const Params p) {
  impl = std::make_shared<Impl>();
  impl->init(mm, m, he, traj, otm, p);
}

void LagrangianLevels
::map_zeta_to_eulerian_cell_centers (const Acr3& p_advect, const Ar3& zeta_vtx,
                                     const Ar1& wrk, const Ar2& zeta_lag) {
  impl->transfer_lv2ec(p_advect, zeta_vtx, wrk, zeta_lag);
}

void LagrangianLevels
::remap_vertically (const Acr2& zeta_eul, const Ar2& zeta_lag, const Ar2& rho,
                    const Ar1& rwrk, const int nq, const Ar3& Q) {
  impl->remap_vertically(zeta_eul, zeta_lag, rho, rwrk, nq, Q);
}

bool LagrangianLevels::unittest (
  const MpasMesh::Ptr& mm, const siqp::sh::Mesh::ConstPtr& m,
  const Comm::Ptr& he, const bool verbose)
{
  bool ok = true;
  // Only 2D flows are used here, so nilvl = nmlvl = nlvl.
  const auto nlvl = mm->l.nmax, ncell = mm->t.ncell, ncown = mm->ncown;
  Ar2 f0_cc(nlvl, ncell), f_cc(nlvl, ncell);
  Ar3 f_cc_3(1, nlvl, ncell, f_cc.data());
  Ar3 plag(3, nlvl, mm->nvown), f_vtx(max_vtx, nlvl, mm->ncown);
  f_vtx.set(-1);
  CellFromCellsModalShapeFunction ccsf;
  ccsf.init(mm, *m, he, Impl::order, Impl::tq_order, false, true);
  Ar1 wrk(OrderInvariantAccumTraits::local_atomic_size_nreal*nlvl*ncell);
  for (const auto vfsym : {'n'}) {
    auto vf = vf_from_short(vfsym);
    set_interface_z_all_cols(*vf, *mm);
    init_zeta_rep(mm, *m, he);
    for (const auto traj_method :
      {Trajectory::Method::exact, Trajectory::Method::taylor2}) {
      const Int nstep = calc_default_nstep(*mm, *he, *vf);
      Int nsubstep = -1, traj_nhalo = -1;
      calc_substep_data(*mm, *he, *vf, nstep, traj_nhalo, nsubstep);
      adjust_nsubstep_nhalo_if_needed(*he, traj_method, 1, nsubstep, traj_nhalo);
      auto traj = std::make_shared<Trajectory>(mm, he, vf, traj_method,
                                               nsubstep, traj_nhalo);
      assert(traj->num_levels() == nlvl);
      Ar1 rwrk(traj->get_rwrk_size());
      OctreeMgr::Ptr otm;
      if ( ! traj->has_in_cell_data())
        otm = std::make_shared<OctreeMgr>(m->p, m->e);
      for (const bool cell_int : {true, false}) {
        LagrangianLevels::Params params;
        params.cell_integrated = cell_int;
        params.bfb = false;
        auto ll = std::make_shared<LagrangianLevels>(mm, m, he, traj, otm, params);
        const Real te = vf->cycle_time(), dt = te/nstep;
        const int tstep = std::max(0, nstep/4 - 1);
        for (int ti = tstep; ti < tstep+1; ++ti) {
          Real l2_prev = 0, ooa;
          for (const int fac : {4, 2}) {
            citest::fill_cc_data(mm, fac, f0_cc);
            const Real ts = (Real(ti)/nstep)*te;
            // Transfer Eul cell data to Eul vtx data.
            ll->impl->transfer_ec2ev(ccsf, f0_cc, f_vtx);
            // Take bwd-in-time step.
            traj->advect_vtxs(m->p, ts, dt, rwrk, plag);
            // Transfer Lag vtx data to Eul cell data.
            ll->impl->transfer_lv2ec(plag, f_vtx, wrk, f_cc);
            // Repeat, this time going forward in time. Need to comm f_cc's halo
            // for the reconstruction in the first step.
            const int ns[] = {1};
            he->exchange_cell_full(1, ns, nlvl, &f_cc_3);
            ll->impl->transfer_ec2ev(ccsf, f_cc, f_vtx);
            traj->advect_vtxs(m->p, ts+dt, -dt, rwrk, plag);
            ll->impl->transfer_lv2ec(plag, f_vtx, wrk, f_cc);
            // Compare final q against original. We expect OOA 2.
            for (Int li = 0; li < 1; ++li) {
              CmpOutput e;
              cmp(*mm, *he, 0, li, nlvl, ncown, nullptr,
                  Ar3(1, nlvl, ncell, f0_cc.data()), nullptr, f_cc_3, e);
              ooa = l2_prev == 0 ? 0 : -std::log(e.l2/l2_prev)/std::log(2);
              if (verbose && he->amroot())
                printf("ut_ll> %d %d %d l2 %1.16e (%1.3e)\n",
                       int(cell_int), li, fac, e.l2, ooa);
              l2_prev = e.l2;
            }
          }
          if (ooa < (cell_int ? 2 : 1.8)) ok = false;
        }
      }
    }
  }
  return ok;
}
