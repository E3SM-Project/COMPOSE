#include "mpas_advect_ir.hpp"
#include "mpas_geometry.hpp"
#include "mpas_ci_util.hpp"
#include "mpas_time_int.hpp"
#include "mpas_siqp.hpp"
#include "mpas_accum.hpp"
#include "mpas_cedr.hpp"
#include "mpas_comm_util.hpp"
#include "mpas_lagrangian_levels.hpp"
#include "mpas_zeta.hpp"

#include <set>
#include <iostream>
#include <sstream>

static int omp_get_thread_num_for_timer () { return mpas_omp_get_thread_num(); }

namespace ir {

struct RemapFormOp::Impl {
  void init(const MpasMesh::ConstPtr& mm, const siqp::sh::Mesh::ConstPtr& m,
            const Comm::Ptr& he, const Params& p);

  void integrate(const Real ts, const Real dt, const Ar2& rhos, const Ar2& rhot,
                 const Ar3& Qs, const Ar3& Qt);

  bool unittest(Real* data = nullptr);
  void advect_mesh(const Real ts, const Real dt);
  void fill_patches(std::vector<Real>& xyzs_orig, std::vector<Real>& xyzs_advected);

  MpasMesh::ConstPtr mm;
  siqp::sh::Mesh::ConstPtr m;
  Comm::Ptr he;
  Params param;

private:
  Int max_vtx, vdegree, nedge, ncown /* nbr own cells */, ntracer;
  // Number of (m)idpoint and (i)nterface levels. nmlvl = nilvl if we're running
  // a set of 2D test problems. nilvl = nmlvl + 1 if the flow is 3D.
  Int nmlvl, nilvl;
  Real avg_diam;
  OctreeMgr::ConstPtr otm;
  Ar3 p_advect; // (3,nilvl,nvown)
  // Trajectory data.
  Trajectory::Ptr trajectory;
  bool have_in_cell_data;
  // For Lagrangian levels.
  bool is_3d;
  LagrangianLevels::Ptr lag_lvls;
  Ar2 dzeta_lag;
  // For high-order reconstruction.
  static constexpr Real rho_lo_thr = 1e-6; // relative to rho ~= 1
  static constexpr int max_ntracer = 15;   // for automatic arrays
  CellFromCellsModalShapeFunction sf;
  Ar4 shapes; // reconstruction data read/written during time stepping
  // For CEDR.
  Ar3 cedr_lo_ini, cedr_hi_ini, cedr_lo_fin, cedr_hi_fin;
  // For upwind MPI.
  bool upwind_mpi;
  //   Does halo cell,lvl (lvl,i) need remote data? Set to >= 1 if in stencil, 2
  //   if potential source.
  Ac2 um_src_needed;
  //   remap_tgt_overlap(0:remap_tgt_noverlap(tc,lvl)-1, lvl, tc) are source
  //   cells that might overlap target cell tc.
  Ai3 remap_tgt_overlap;
  Ai2 remap_tgt_noverlap;
  // Workspace used by multiple steps that don't overlap data usage.
  Ar1 rbuf;

  // (lvl,sc) index space.
  template<typename Fn> inline void ompparfor_cell_lvl(const Fn& f, const Int ncell_max);
  template<typename Fn> inline void ompfor_cell_lvl(const Fn& f, const Int ncell_max);

  // Exchange source cell data.
  void exchange_src_tracer_data(const Ar2& rhos, const Ar3& Qs);

  // If the flow is 3D, handle vertical part of time integration.
  void remap_src_vertically(const Ar2& rhos, const Ar3& Qs);

  // Reconstruct the source data.
  void reconstruct(const Acr2& rhos, const Acr3& Qs);

  // Remap from source to target.
  void remap(const TriangleQuadrature& tq, const Ar2& rhos, const Ar2& rhot,
             const Ar3& Qs, const Ar3& Qt);
  void search(const Int tid, const Int tc, const Int lvl, OTSearchFunctor& f);
  void create_advected_polygon(const Int tc, const Int lvl, PolyClipper& cd);
  void accumulate_integrals(const int tid, const TriangleQuadrature& tq,
                            const Acr2& rhos, const Acr3& Qs, const Int sc,
                            const Int lvl, const PolyClipper& pc, Real* accum);
  void update_bounds(const Int ncell, const Int sc, const Int tc, const Int lvl);
  void set_tgt_values(const Int ncell, const Int tc, const Int lvl, const Real* accum,
                      const Ar2& rhot, const Ar3& Qt);
  void run_cedr(const Acr2& rhos, const Ar2& rhot, const Acr3& Qs, const Ar3& Qt);
};

template<typename Fn>
void RemapFormOp::Impl::ompfor_cell_lvl (const Fn& f, const Int ncell_max) {
  const Int n = ncell_max*nmlvl;
  const auto& l = mm->l;
  ompfor for (Int i = 0; i < n; ++i) {
    const Int ci = i / nmlvl, li = i % nmlvl;
    if (li < l.get_cmin(ci) || li > l.get_cmax(ci)) continue;
    f(ci, li);
  }
}

template<typename Fn>
void RemapFormOp::Impl::ompparfor_cell_lvl (const Fn& f, const Int ncell_max) {
  omppar ompfor_cell_lvl(f, ncell_max);
}

void RemapFormOp::Impl
::init (const MpasMesh::ConstPtr& mm_, const siqp::sh::Mesh::ConstPtr& m_,
        const Comm::Ptr& he_, const Params& p_) {
  mm = mm_;
  m = m_;
  he = he_;
  param = p_;
  if (param.tq_order <= 0) param.tq_order = param.cedr ? 2 : 8;
  ncown = mm->ncown;
  assert(ncown == get_num_own_cells(mm->t, *he));
  max_vtx = mm->t.max_vtx;
  vdegree = mm->t.vdegree;
  nedge = mm->t.nedge;
  ntracer = param.ntracer;
  is_3d = param.vf->is_3d();
  nmlvl = mm->l.nmax;
  nilvl = nmlvl + (is_3d ? 1 : 0);
  upwind_mpi = he->nrank() > 1;

  p_advect.reset(is_3d ? 4 : 3, nilvl, mm->nvown);
  size_t rwrk = 0;

  sf.init(mm, *m, he, param.order, param.tq_order, param.limit,
          false /* owned and shared */);
  if (param.order >= 2) {
    // Shapes data for fields.
    rwrk = size_t(sf.ho_ncoef)*(1 + ntracer)*nmlvl*mm->t.ncell;
  }

  adjust_nsubstep_nhalo_if_needed(*he, param.advect_method,
                                  CellShapeFunction::two_halo ? 2 : 1,
                                  param.nsubstep, param.traj_nhalo);
  trajectory = std::make_shared<Trajectory>(mm, he, param.vf, param.advect_method,
                                            param.nsubstep, param.traj_nhalo);
  have_in_cell_data = (trajectory->has_in_cell_data() && ! param.perturb_dep_pts);
  if ( ! have_in_cell_data)
    otm = std::make_shared<OctreeMgr>(m->p, m->e);
  rwrk = std::max(rwrk, trajectory->get_rwrk_size());

  if (is_3d) {
    LagrangianLevels::Params p;
    p.cell_integrated = param.lag_lvls_cell_int;
    p.bfb = param.bfb;
    lag_lvls = std::make_shared<LagrangianLevels>(mm, m, he, trajectory, otm, p);
    dzeta_lag.reset(nilvl, mm->t.ncell);
    // We share some memory with shapes.
    rwrk = std::max(rwrk, size_t(nilvl)*max_vtx*ncown);
    // For Trajectory::remap_vertically.
    rwrk = std::max(rwrk, size_t(nmlvl)*4*mpas_omp_get_max_threads());
  }

  if (param.cedr)
    rwrk = std::max<size_t>(
      rwrk, (param.bfb ?
             // See caas docs for details.
             (OrderInvariantAccumTraits::global_size_nreal*ntracer*
              (ncown + mpas_omp_get_max_threads() + 1)) :
             6*ntracer));

  if (is_3d && param.lag_lvls_cell_int && param.bfb)
    rwrk = std::max<size_t>(
      rwrk, (OrderInvariantAccumTraits::local_atomic_size_nreal*
             nilvl*mm->t.ncell));

  if (rwrk) {
    rbuf.reset(rwrk);
    if (param.order >= 2)
      shapes.reset(sf.ho_ncoef, 1 + ntracer, nmlvl, mm->t.ncell, rbuf.data());
  }

  Real radius;
  get_avg_radius_and_cell_spacing(*mm, *he, radius, avg_diam);
}

void RemapFormOp::Impl
::reconstruct (const Acr2& rhos, const Acr3& Qs) {
  const int ntracer = param.ntracer;
  const int ncell = m->e.n();
  const auto ho_ncoef = sf.ho_ncoef;
  Acr3 rhos3(1, nmlvl, ncell, rhos.data());
  const auto f = [&] (const Int sc, const Int li) {
    if (at_halo_edge(*he, sc)) return;
    if (upwind_mpi && um_src_needed(li,sc) < 2) return;
    sf.calc_shape(sc, rhos3, 0, li, shapes(0,li,sc));
    for (int qi = 0; qi < ntracer; ++qi)
      sf.calc_shape(sc, Qs, qi, li, shapes(qi+1,li,sc));
  };
  ompparfor_cell_lvl(f, ncell);
  if (param.limit || param.cedr) {
    const auto g = [&] (const Int sc, const Int li) {
      if (at_halo_edge(*he, sc)) return;
      if (upwind_mpi && um_src_needed(li,sc) < 2) return;
      const int* const stencil = sf.c2c(sc);
      const int stensz = sf.cell_nstencil[sc];
      for (int qi = 0; qi < ntracer; ++qi) {
        Real qs_min = Qs(qi,li,sc)/rhos(li,sc), qs_max = qs_min;
        for (int i = 0; i < stensz; ++i) {
          const int cj = stencil[i];
          const auto q = Qs(qi,li,cj)/rhos(li,cj);
          qs_min = std::min(qs_min, q);
          qs_max = std::max(qs_max, q);
        }
        if (param.cedr) {
          cedr_lo_ini(qi,li,sc) = qs_min;
          cedr_hi_ini(qi,li,sc) = qs_max;
        }
        if ( ! param.limit) continue;
        assert(param.order == 2); // limiter works only for linear shape function
        assert( ! CellShapeFunction::scale_by_arc_length);
        const Real* const cc = sf.ho_coefs(sc);
        const Real* const xaxis = cc + 3;
        const Real* const yaxis = xaxis + 3;
        const int ne = mm->t.c2ne[sc];
        assert(ho_ncoef == 2);
        Real rho_shape[2];
        geo::copyn(ho_ncoef, rho_shape, shapes(0, li, sc));
        if (param.cedr && qi == 0)
          for (int i = 0; i < ne; ++i) {
            const auto v = m->p(m->e(i,sc));
            Vec vmcc; geo::axpbyz(1, v, -1, cc, vmcc);
            const Real x = geo::dot(vmcc, xaxis), y = geo::dot(vmcc, yaxis);
            const Real ho_term = rho_shape[0]*x + rho_shape[1]*y;
            const Real rhos_reconstruct = rhos(li,sc) + ho_term;
            if (rhos_reconstruct < rho_lo_thr) {
              Real alpha = 0;
              if (rhos(li,sc) >= rho_lo_thr)
                alpha = (rho_lo_thr - rhos(li,sc))/ho_term;
              rho_shape[0] *= alpha; rho_shape[1] *= alpha;
            }
          }
        Real shape[2];
        geo::copyn(ho_ncoef, shape, shapes(qi+1, li, sc));
        for (int i = 0; i < ne; ++i) {
          const auto v = m->p(m->e(i,sc));
          Vec vmcc; geo::axpbyz(1, v, -1, cc, vmcc);
          const Real
            x = geo::dot(vmcc, xaxis),
            y = geo::dot(vmcc, yaxis),
            rho_lo = rhos(li,sc),
            rho_hi = rho_lo + rho_shape[0]*x + rho_shape[1]*y,
            Q_lo = Qs(qi,li,sc),
            Q_hi = Q_lo + shape[0]*x + shape[1]*y;
          const bool below = Q_hi < qs_min*rho_hi, above = Q_hi > qs_max*rho_hi;
          if (below || above) {
            // Q_reconstruct = alpha Q_hi + (1 - alpha) (Q_lo/rho_lo) rho_hi,
            // where alpha is chosen to put Q_reconstruct at the bounds.
            const Real q_lo = Q_lo/rho_lo, q_hi = Q_hi/rho_hi;
            Real alpha = (q_lo == q_hi ?
                          0 :
                          (below ?
                           (qs_min - q_lo)/(q_hi - q_lo):
                           (qs_max - q_lo)/(q_hi - q_lo)));
            if (alpha < 0) alpha = 0;
            for (int j = 0; j < 2; ++j)
              shape[j] = alpha*shape[j] + (1 - alpha)*q_lo*rho_shape[j];
          }
        }
        if (qi == 0) geo::copyn(ho_ncoef, shapes(0, li, sc), rho_shape);
        geo::copyn(ho_ncoef, shapes(qi+1, li, sc), shape);
      }
    };
    ompparfor_cell_lvl(g, ncell);
  }
}

void RemapFormOp::Impl
::search (const Int tid, const Int tc, const Int lvl, OTSearchFunctor& f) {
  if (have_in_cell_data) {
    const auto& t = mm->t;
    const int ne = t.c2ne[tc];
    assert(f.nhit == 0 && f.ok);
    for (int i = 0; i < ne; ++i) {
      int cis[2];
      cis[0] = trajectory->
        get_cell_lvl_containing_vtx(t.get_c2v(tc,i), lvl).cell();
      if (is_3d)
        cis[1] = trajectory->
          get_cell_lvl_containing_vtx(t.get_c2v(tc,i), lvl+1).cell();
      for (int i = 0; i < (is_3d ? 2 : 1); ++i) {
        const int ci = cis[i];
        if (i == 1 && ci == cis[0]) break;
        f(cis[i]);
        if (at_halo_edge(*he, ci)) {
          // Too close to the halo edge. This will be caught later. For now, stop
          // searching.
          continue;
        }
        assert(ci >= 0 && ci < t.ncell);
        const int ine = t.c2ne[ci];
        for (int j = 0; j < ine; ++j) {
          const auto cj = t.get_c2c(ci,j);
          assert(cj >= 0 && cj < t.ncell);
          f(cj);
        }
      }
    }
    assert(f.ok); // should be plenty of room in HitsAccumulator
  } else {
    Timer::start(tid, Timer::search);
    if (is_3d) {
      const int ne = mm->t.c2ne[tc];
      int e[::max_vtx];
      for (int i = 0; i < ne; ++i) e[i] = i;
      Real buf[3*::max_vtx];
      const Ar3 p(3, 1, max_vtx, buf);
      const auto e_advect = m->e(tc);
      for (int i = 0; i < ne; ++i)
        for (int d = 0; d < 3; ++d)
          p(d,0,i) = 0.5*(p_advect(d,lvl  ,e_advect[i]) +
                          p_advect(d,lvl+1,e_advect[i]));
      otm->apply(const_unmgd(p), 0, e, ne, f);
    } else {
      otm->apply(const_unmgd(p_advect), lvl, m->e(tc), m->e.m(), f);
    }
    // If this happens, something else has already gone wrong in the numerics.
    require_msg(f.ok, "search !have_in_cell_data: Too many hits.");
    Timer::stop(tid, Timer::search);
  }
  const Int* const gids = get_cell_gids(*he);
  f.post(gids);
}

void RemapFormOp::Impl
::create_advected_polygon (const Int tc, const Int lvl, PolyClipper& cd) {
  const int ne = mm->t.c2ne[tc];
  cd.ni = ne;
  if (nilvl == nmlvl) {
    // We're running a 2D-flow test problem, but still possibly with multiple
    // levels.
    const auto e = m->e(tc);
    for (int i = 0; i < ne; ++i)
      geo::copy(cd.vi(i), p_advect(lvl, e[i]));
  } else {
    // We're running a 3D-flow problem.
    const auto e = m->e(tc);
    for (int i = 0; i < ne; ++i)
      geo::axpbyz(0.5, p_advect(lvl  ,e[i]),
                  0.5, p_advect(lvl+1,e[i]), cd.vi(i));
  }
}

void RemapFormOp::Impl
::accumulate_integrals (const int tid, const TriangleQuadrature& tq,
                        const Acr2& rhos, const Acr3& Qs, const Int sc,
                        const Int lvl, const PolyClipper& pc, Real* accum) {
  Timer::start(tid, Timer::quadrature);
  Real basis[CellShapeFunction::max_ho_ncoef+1];
  sf.integrate_basis_over_polygon(tq, pc.vo, pc.no, sc, basis);
  Real dz = 1;
  if (is_3d) dz = dzeta_lag(lvl,sc);
  if (param.order == 1) {
    accum[0] += sf.integrate(rhos(lvl,sc), nullptr, basis)*dz;
    for (int qi = 0; qi < ntracer; ++qi)
      accum[qi+1] += sf.integrate(Qs(qi,lvl,sc), nullptr, basis)*dz;
  } else {
    accum[0] += sf.integrate(rhos(lvl,sc), shapes(0,lvl,sc), basis)*dz;
    for (int qi = 0; qi < ntracer; ++qi)
      accum[qi+1] += sf.integrate(Qs(qi,lvl,sc), shapes(qi+1,lvl,sc), basis)*dz;
  }
  Timer::stop(tid, Timer::quadrature);
}

void RemapFormOp::Impl
::update_bounds (const Int ncell, const Int sc, const Int tc, const Int lvl) {
  for (int qi = 0; qi < ntracer; ++qi) {
    // Use the full halo of the reconstruction.
    cedr_lo_fin(qi,lvl,tc) = std::min(cedr_lo_fin(qi,lvl,tc), cedr_lo_ini(qi,lvl,sc));
    cedr_hi_fin(qi,lvl,tc) = std::max(cedr_hi_fin(qi,lvl,tc), cedr_hi_ini(qi,lvl,sc));
  }
}

void RemapFormOp::Impl
::set_tgt_values (const Int ncell, const Int tc, const Int lvl, const Real* accum,
                  const Ar2& rhot, const Ar3& Qt) {
  Real volume = mm->g.my_areaCell[tc];
  if (is_3d) {
    const auto& zeta_eul = mm->zeta_rep->eul_cell;
    volume *= zeta_eul(lvl+1,tc) - zeta_eul(lvl,tc);
  }
  rhot(lvl,tc) = accum[0]/volume;
  for (int qi = 0; qi < ntracer; ++qi) Qt(qi,lvl,tc) = accum[qi+1]/volume;
}

static inline void check_sc (const Comm& he, const Int sc) {
  require_msg( ! at_halo_edge(he, sc),
              "Source cell is on comm halo edge; "
              "increase comm halo or decrease time step.");  
}

void RemapFormOp::Impl
::remap (const TriangleQuadrature& tq, const Ar2& rhos, const Ar2& rhot,
         const Ar3& Qs, const Ar3& Qt) {
  Timer::start(Timer::scq);
  const int ncell = m->e.n();

  const auto f = [&] (const Int tc, const Int li) {
    const int tid = omp_get_thread_num_for_timer();

    assert(ntracer <= max_ntracer);
    Real accum[max_ntracer+1] = {0};

    if (upwind_mpi) {
      PolyClipper pc;
      create_advected_polygon(tc, li, pc);
      const Int noverlap = remap_tgt_noverlap(li,tc);
      const Int* const overlap = remap_tgt_overlap(li,tc);
      for (Int i = 0; i < noverlap; ++i) {
        const Int sc = overlap[i];
        assert(sc >= 0 && sc < ncell);
        check_sc(*he, sc);
        pc.clip(*m, sc);
        assert(pc.no);
        accumulate_integrals(tid, tq, rhos, Qs, sc, li, pc, accum);
        if (param.cedr) update_bounds(ncell, sc, tc, li);
      }
    } else {
      // Find potential overlapping source cells.
      OTSearchFunctor f;
      search(tid, tc, li, f);
      // Clip target against potential source cell and, if non-empty, integrate
      // over the overlap.
      PolyClipper pc;
      create_advected_polygon(tc, li, pc);
      for (int hi = 0; hi < f.nhit; ++hi) {
        const int sc = f.hits[hi];
        if (sc == -1) continue;
        check_sc(*he, sc);
        pc.clip(*m, sc);
        if (pc.no) { // If no > 0, an overlap was found, so accumulate.
          accumulate_integrals(tid, tq, rhos, Qs, sc, li, pc, accum);
          if (param.cedr) update_bounds(ncell, sc, tc, li);
        }
      }
      // If this happens, something else has already gone wrong in the numerics.
      require_msg(f.ok, "remap !upwind_mpi: Too many hits.");
    }

    set_tgt_values(ncell, tc, li, accum, rhot, Qt);
  };
  ompparfor_cell_lvl(f, ncown);
  Timer::stop(Timer::scq);
}

void RemapFormOp::Impl
::run_cedr (const Acr2& rhos, const Ar2& rhot, const Acr3& Qs, const Ar3& Qt) {
  assert(param.cedr);
  Timer::start(Timer::cedr);

  const auto& zeta_eul = mm->zeta_rep->eul_cell;

  // Make target rho positive and conserve mass in 3D case. In practice rho is
  // obtained externally, in which case this all_reduce is avoided.
  if (is_3d) {
    const auto f = [&] (const int ci, const int li) {
      if (rhot(li,ci) < rho_lo_thr) rhot(li,ci) = rho_lo_thr;
    };
    ompparfor_cell_lvl(f, ncown);
    const auto g = [&] (int ci, Real a[2]) {
      // With Intel 19.1.3.304 and -fp-model fast (default), 'volatile' is
      // needed in the decls of a0 and a1 to get PE-invariant
      // BFBness. Alternatively, with -fp-model source, precise, etc applied to
      // just run_cedr in its own cpp file, BFBness is obtained.
      Real volatile a0 = 0, a1 = 0;
      for (Int li = 0; li < nmlvl; ++li) {
        const Real dzs = dzeta_lag(li,ci);
        const Real dzt = zeta_eul(li+1,ci) - zeta_eul(li,ci);
        a0 += rhos(li,ci)*dzs;
        a1 += rhot(li,ci)*dzt;
      }
      const Real area = mm->g.my_areaCell[ci];
      a[0] = area*a0; a[1] = area*a1;
    };
    Real gm[2];
    repro_sum(*he, ncown, 2, rbuf, g, param.bfb, gm);
    const Real fac = gm[0]/gm[1];
    const auto h = [&] (const int ci, const int li) { rhot(li,ci) *= fac; };
    ompparfor_cell_lvl(h, ncown);
  } else {
    for (Int li = 0; li < nmlvl; ++li) {
      const auto f = [&] (int ci, Real a[2]) {
        a[0] = rhos(li,ci)*mm->g.my_areaCell[ci];
        if (rhot(li,ci) < rho_lo_thr) rhot(li,ci) = rho_lo_thr;
        a[1] = rhot(li,ci)*mm->g.my_areaCell[ci];
      };
      Real gm[2];
      repro_sum(*he, ncown, 2, rbuf, f, param.bfb, gm);
      const Real fac = gm[0]/gm[1];
      ompparfor for (int ci = 0; ci < ncown; ++ci) rhot(li,ci) *= fac;
    }
  }

  { // Bounds on target Q.
    const auto f = [&] (const int ci, const int li) {
      const auto rho = rhot(li,ci);
      for (int qi = 0; qi < ntracer; ++qi) {
        cedr_lo_fin(qi,li,ci) *= rho;
        cedr_hi_fin(qi,li,ci) *= rho;
      }
    };
    ompparfor_cell_lvl(f, ncown);
  }

  // Run CAAS. N.B. that a single problem includes all levels.
  caas(*he, ntracer, nmlvl, ncown, mm->g.my_areaCell.data(), is_3d, dzeta_lag,
       zeta_eul, cedr_lo_fin, cedr_hi_fin, Qs, Qt, Qt, rbuf, false,
       param.bfb);

  Timer::stop(Timer::cedr);
}

void RemapFormOp::Impl
::advect_mesh (const Real ts, const Real dt) {
  Timer::start(Timer::trajectory);
  trajectory->advect_vtxs(m->p, ts, dt, rbuf, p_advect);
  Timer::stop(Timer::trajectory);
  if (param.perturb_dep_pts) {
    // Perturb the departure points to show that with if we apply a limiter in
    // the source reconstruction and a CEDR to the target, we are robust to cell
    // self-intersections.
    const auto perturb = [&] (Real* const p_a) {
      for (int d = 0; d < 3; ++d) p_a[d] += avg_diam*(urand() - 0.5);
      geo::normalize(p_a);
    };
    for (Int li = 0; li < nilvl; ++li)
      for (Int i = 0; i < mm->nvown; ++i)
        perturb(p_advect(li,i));
  }
}

void RemapFormOp::Impl
::remap_src_vertically (const Ar2& rhos, const Ar3& Qs) {
  if ( ! is_3d) return;
  // In this routine, dzeta_lag is actually zeta_lag until the last step.
  const auto ncown = mm->ncown;
  // Remap zeta to Eulerian cell centers.
  {
    const Ar3 zeta_vtx(max_vtx, nilvl, ncown, rbuf.data());
    const auto& t = mm->t;
    assert(mm->zeta_rep);
    const auto& eul_vtx_multivalued = mm->zeta_rep->eul_vtx_multivalued;
    const auto& eul_vtx_univalued = mm->zeta_rep->eul_vtx_univalued;
    ompparfor for (int ci = 0; ci < ncown; ++ci) {
      const int ne = t.c2ne[ci];
      const auto e = m->e(ci);
      for (int li = 0; li < nilvl; ++li)
        for (int i = 0; i < ne; ++i) {
          zeta_vtx(i,li,ci) =
            // Combine the good but vertex-multivalued approximation to the
            // Eulerian zeta with the increment in the vertex-univalued but not
            // as good approximation.
            //   In the case that zeta_eul is constant within a level, this
            // arithmetic simplifies to just using p_advect(3,:,:) directly.
            //   Crucially, as dt->0, this approach assures that zeta_vtx ->
            // zeta_eul; thus, there is no O(1) error due simply to remap.
            ( // Good approximation to possibly spatially variable Eulerian
              // zeta, zeta_eul.
              eul_vtx_multivalued(i,li,ci) +
              // Delta zeta for univalued approximation of zeta_eul, which is
              // not as good as the multivalued one.
              (p_advect(3,li,e[i]) - eul_vtx_univalued(li,e[i])));
        }
    }
    lag_lvls->map_zeta_to_eulerian_cell_centers(p_advect, zeta_vtx, rbuf,
                                                dzeta_lag);
  }
  // Vertically Eulerian grid.
  const auto& zeta_eul = mm->zeta_rep->eul_cell;
  // Explicitly enforce exactness of src and tgt domains.
  ompparfor for (int ci = 0; ci < ncown; ++ci) {
    dzeta_lag(0,ci) = zeta_eul(0,ci);
    dzeta_lag(nilvl-1,ci) = zeta_eul(nilvl-1,ci);
  }
  // Vertically remap rho and Q at cell centers.
  lag_lvls->remap_vertically(zeta_eul, dzeta_lag, rhos, rbuf, ntracer, Qs);
  // Diff zeta_lag in the vertical.
  ompparfor for (int ci = 0; ci < ncown; ++ci)
    for (int li = 0; li < nmlvl; ++li)
      dzeta_lag(li,ci) = dzeta_lag(li+1,ci) - dzeta_lag(li,ci);
}

void RemapFormOp::Impl
::exchange_src_tracer_data (const Ar2& rhos, const Ar3& Qs) {
  const Int ncell = m->e.n();
  if (upwind_mpi) {
    if (remap_tgt_overlap.empty()) {
      remap_tgt_overlap.reset(PolyClipper::max_src_overlap, nmlvl, ncown);
      remap_tgt_noverlap.reset(nmlvl, ncown);
    }
    if (um_src_needed.empty()) um_src_needed.reset(nmlvl,ncell);
    ompparfor for (Int i = 0; i < ncell; ++i)
      for (Int j = 0; j < nmlvl; ++j)
        um_src_needed(j,i) = 0;
    const auto f = [&] (const Int tc, const Int li) {
      const int tid = omp_get_thread_num_for_timer();
      OTSearchFunctor f;
      search(tid, tc, li, f);
      PolyClipper pc;
      create_advected_polygon(tc, li, pc);
      remap_tgt_noverlap(li,tc) = 0;
      for (int hi = 0; hi < f.nhit; ++hi) {
        const int sc = f.hits[hi];
        if (sc == -1) continue;
        check_sc(*he, sc);
        pc.clip(*m, sc);
        if (pc.no) { // If no > 0, an overlap was found, so record this src cell.
          um_src_needed(li,sc) = 2;
          // Record overlap for use in remap.
          remap_tgt_overlap(remap_tgt_noverlap(li,tc),li,tc) = sc;
          ++remap_tgt_noverlap(li,tc);
          assert(remap_tgt_noverlap(li,tc) <= pc.max_src_overlap);
        }
      }
    };
    ompparfor_cell_lvl(f, ncown);
    if (param.order > 1) {
      // Need data for nbr cells of cells needing reconstruction.
      const auto g = [&] (const Int sc, const Int li) {
        if (um_src_needed(li,sc) != 2) return;
        const int stensz = sf.cell_nstencil[sc];
        const int* stencil = sf.c2c(sc);
        for (Int j = 0; j < stensz; ++j) {
          if (stencil[j] < ncown) continue;
          assert(stencil[j] < ncell);
          auto& usn = um_src_needed(li,stencil[j]);
          // The next lines can convert 0 to 1 and 1 to 1 but can't change 2.
          if (usn == 2) continue;
          usn = 1; // benign write race
        }
      };
      ompparfor_cell_lvl(g, ncell);
    }
    const Ar3 rhos3(1, nmlvl, ncell, const_cast<Real*>(rhos.data()));
    const char maskval = 1 | 2;
    if (is_3d) {
      const Ar3 dzeta_lag3(1, nilvl, ncell, const_cast<Real*>(dzeta_lag.data()));
      int ns[] = {1, 1, ntracer};
      const Ar3 fs[] = {dzeta_lag3, rhos3, Qs};
      he->exchange_cell_masked(
        Array2D<const char>(nmlvl, ncell-ncown, um_src_needed(ncown)),
        maskval, 3, ns, nmlvl, fs);
    } else {
      int ns[] = {1, ntracer};
      const Ar3 fs[] = {rhos3, Qs};
      he->exchange_cell_masked(
        Array2D<const char>(nmlvl, ncell-ncown, um_src_needed(ncown)),
        maskval, 2, ns, nmlvl, fs);
    }
  } else {
    const Ar3 rhos3(1, nmlvl, ncell, const_cast<Real*>(rhos.data()));
    if (is_3d) {
      const Ar3 dzeta_lag3(1, nilvl, ncell, const_cast<Real*>(dzeta_lag.data()));
      int ns[] = {1, 1, ntracer};
      const Ar3 fs[] = {dzeta_lag3, rhos3, Qs};
      he->exchange_cell_full(3, ns, nmlvl, fs);
    } else {
      int ns[] = {1, ntracer};
      const Ar3 fs[] = {rhos3, Qs};
      he->exchange_cell_full(2, ns, nmlvl, fs);
    }
  }
}

void RemapFormOp::Impl
::integrate (const Real ts, const Real dt, const Ar2& rhos, const Ar2& rhot,
             const Ar3& Qs, const Ar3& Qt) {
  Timer::start(Timer::integrate);

  const int ntracer = param.ntracer;
  const int ncell = mm->t.ncell;

  if (param.cedr) {
    cedr_lo_ini.reset(ntracer,nmlvl,ncell); cedr_hi_ini.reset(ntracer,nmlvl,ncell);
    cedr_lo_fin.reset(ntracer,nmlvl,ncown); cedr_hi_fin.reset(ntracer,nmlvl,ncown);
  }

  advect_mesh(ts, dt);

  remap_src_vertically(rhos, Qs);

  exchange_src_tracer_data(rhos, Qs);

  TriangleQuadrature tq;
  Timer::start(Timer::reconstruct);
  if (param.order > 1) reconstruct(rhos, Qs);
  Timer::stop(Timer::reconstruct);

  if (param.cedr) {
    omppar {
      if (param.order == 1) {
        const auto g = [&] (const int ci, const int li) {
          for (int qi = 0; qi < ntracer; ++qi) {
            if (at_halo_edge(*he, ci)) continue;
            if (upwind_mpi && ! um_src_needed(li,ci)) continue;
            cedr_lo_ini(qi,li,ci) = cedr_hi_ini(qi,li,ci) = Qs(qi,li,ci)/rhos(li,ci);
          }
        };
        ompfor_cell_lvl(g, ncell);
      }
      const auto f = [&] (const int ci, const int li) {
        for (int qi = 0; qi < ntracer; ++qi) {
          cedr_lo_fin(qi,li,ci) =  1e30;
          cedr_hi_fin(qi,li,ci) = -1e30;
        }
      };
      ompfor_cell_lvl(f, ncown);
    }
  }

  remap(tq, rhos, rhot, Qs, Qt);

  if (param.cedr) run_cedr(rhos, rhot, Qs, Qt);

  Timer::stop(Timer::integrate);
}

RemapFormOp::Ptr RemapFormOp
::create (const MpasMesh::ConstPtr& mm, const siqp::sh::Mesh::ConstPtr& m,
          const Comm::Ptr& he, const Params& p) {
  return std::make_shared<RemapFormOp>(mm, m, he, p);
}

bool RemapFormOp::Impl::unittest (Real* data) {
  return trajectory->unittest(m->p, data) && Trajectory::convtest(mm, he);
}

RemapFormOp
::RemapFormOp (const MpasMesh::ConstPtr& mm, const siqp::sh::Mesh::ConstPtr& m,
               const Comm::Ptr& he, const Params& p) {
  impl = std::make_shared<Impl>();
  impl->init(mm, m, he, p);
}

// If is_3d, rhos and Qs are modified on output.
void RemapFormOp
::integrate (const Real ts, const Real dt, const Ar2& rhos, const Ar2& rhot,
             const Ar3& Qs, const Ar3& Qt)
{ impl->integrate(ts, dt, rhos, rhot, Qs, Qt); }
bool RemapFormOp::unittest (Real* data) { return impl->unittest(data); }
void RemapFormOp::advect_mesh (const Real ts, const Real dt)
{ impl->advect_mesh(ts, dt); }
void RemapFormOp
::fill_patches (std::vector<Real>& xyzs_orig, std::vector<Real>& xyzs_advected)
{ impl->fill_patches(xyzs_orig, xyzs_advected); }

MpasMesh::ConstPtr RemapFormOp::get_mpas_mesh () const { return impl->mm; }
siqp::sh::Mesh::ConstPtr RemapFormOp::get_siqp_mesh () const { return impl->m; }
Comm::Ptr RemapFormOp::get_comm () const { return impl->he; }
const Params& RemapFormOp::get_params () const { return impl->param; }

void RemapFormOp::Impl
::fill_patches (std::vector<Real>& xyzs_orig, std::vector<Real>& xyzs_advected) {
  const int
    ncell = mm->t.ncell,
    nper1d = mm->t.max_vtx, nper3d = 3*nper1d;
  xyzs_orig.resize(ncell*nper3d);
  xyzs_advected.resize(ncell*nper3d);
  ompparfor for (int ci = 0; ci < ncell; ++ci) {
    Real
      * xso = xyzs_orig.data() + ci*nper3d,
      * yso = xso + nper1d,
      * zso = yso + nper1d,
      * xsa = xyzs_advected.data() + ci*nper3d,
      * ysa = xsa + nper1d,
      * zsa = ysa + nper1d;
    const int ne = mm->t.c2ne[ci];
    const auto elem = m->e(ci);
    for (int i = 0; i < ne; ++i) {
      auto p = m->p(elem[i]);
      xso[i] = p[0]; yso[i] = p[1]; zso[i] = p[2];
      p = p_advect(0, elem[i]);
      xsa[i] = p[0]; ysa[i] = p[1]; zsa[i] = p[2];
    }
  }
}

} // namespace ir
