#include "mpas_trajectory_impl.hpp"

#include <set>

typedef Trajectory::CellLvl CellLvl;

// Iterate over entities of type cell, edge, or vtx and levels.
template<typename Fn>
static void ompfor_ent_lvl (const Fn& f, const Levels& l, const bool is_3d,
                            const Int n_max) {
  const Int delta = is_3d ? 1 : 0;
  const Int nlvl = l.nmax + delta, n = n_max*nlvl;
  ompfor for (Int i = 0; i < n; ++i) {
    const Int ci = i / nlvl, li = i % nlvl;
    if (li < l.get_cmin(ci) || li > l.get_cmax(ci) + delta) continue;
    f(ci, li);
  }
}

template<typename Fn>
static void ompparfor_ent_lvl (const Fn& f, const Levels& l, const bool is_3d,
                               const Int n_max) {
  omppar ompfor_ent_lvl(f, l, is_3d, n_max);
}

enum struct MaskVals : char {
  from_cell     = 1 << 0,
  from_cell_it1 = 1 << 1,
  from_vtx      = 1 << 6
};

static void mask_zero (const Ac2& mask) {
  mask.set_ompfor(0);
}

static void mask_and (const Ac2& mask, const char value) {
  const Int m = mask.m(), n = mask.n();
  ompfor for (int i = 0; i < n; ++i)
    for (int j = 0; j < m; ++j)
      mask(j,i) &= value;
}

Trajectory::Method::Enum
Trajectory::Method::from_short (const char sym) {
  switch (sym) {
  case 'e': return exact;
  case '1': return order1;
  case '2': return order2;
  case 't': return taylor2;
  default:
    printf("Method: %c is not valid\n", sym);
    assert(0);
    return exact;
  }
}

const char* Trajectory::Method::to_string (const Enum e) {
  switch (e) {
  case exact:   return "exact";
  case order1:  return "uv-edge-nml-w-cc-1";
  case order2:  return "uv-edge-nml-w-cc-2";
  case taylor2: return "uv-edge-nml-w-cc-taylor-series-2";
  default: assert(0); return "none";
  }
}

Trajectory
::Trajectory (const MpasMesh::ConstPtr& mm, const Comm::Ptr& he, VelocityField::Ptr vf,
              const Method::Enum method, const int nsubstep, const int nhalo) {
  impl = std::make_shared<Impl>();
  impl->init(mm, he, vf, method, nsubstep, nhalo);
}

void Trajectory::set_method_exact_rel_tol (const Real reltol) {
  impl->method_exact_rel_tol = reltol;
}

int Trajectory::num_levels () const { return impl->nlvl; }

void Trajectory::Impl
::init (const MpasMesh::ConstPtr& mm_, const Comm::Ptr& he_, VelocityField::Ptr vf_,
        const Method::Enum method_, const int nsubstep_, const int nhalo_) {
  method = method_;
  nsubstep = nsubstep_;
  nhalo = Method::is_exact(method) ? 1 : nhalo_;
  vf = vf_;
  mm = mm_;
  he = he_;
  is_3d = vf->is_3d();
  nlvl = mm->l.nmax + (is_3d ? 1 : 0);
  method_exact_rel_tol = 1e8*std::numeric_limits<Real>::epsilon();
  require_msg( ! he->parallel() || Method::is_exact(method) ||
              nsubstep*nhalo >= he->szhalo() - 1,
              "Trajectory: if MPI is used, then nsubstep*nhalo must be >= "
               "he->szhalo()-1;" << " nsubstep " << nsubstep << " nhalo "
               << nhalo << " he->szhalo() " << he->szhalo());
  vdegree = mm->t.vdegree;
  setup_advect_vtxs();
  if ( ! Method::is_exact(method)) {
    const auto& t = mm->t;
    vtx_mask.reset(nlvl, mm->t.nvtx);
    edge_mask.reset(nlvl, mm->t.nedge - mm->neown);
    cell_mask.reset(nlvl, mm->t.ncell - mm->ncown);
    // Set the mask for substep 0, iteration 0, at init time, and never zero
    // that bit.
    const bool cell_halo = method == Method::taylor2;
    omppar {
      mask_zero(edge_mask);
      mask_zero(vtx_mask);
      mask_zero(cell_mask);
      masks_init_from_own((char) MaskVals::from_vtx, cell_halo);
    }
  }
}

void Trajectory::Impl
::masks_init_from_own (const char maskval, const bool cell_halo) {
  const auto& t = mm->t;
  const auto& l = mm->l;
  const auto nvown = mm->nvown, neown = mm->neown, ncown = mm->ncown;
  const auto nedge = t.nedge;
  const auto e2v_nedge = e2v_stencils.m();
  // All of the mask-write lines are benign write race conditions.
  const auto f = [&] (const int vi, const int li) {
    vtx_mask(li,vi) |= maskval;
    const auto* const stencil = e2v_stencils(vi);
    for (Int j = 0; j < e2v_nedge; ++j) {
      const Int ej = stencil[j];
      if (ej >= neown) {
        assert(ej < nedge);
        edge_mask(li,ej-neown) |= maskval;
      }
    }
    if (is_3d) {
      for (Int j = 0; j < t.vdegree; ++j) {
        const Int cj = t.get_v2c(vi,j);
        if (cj >= ncown)
          cell_mask(li,cj-ncown) |= maskval;
      }
    }
    if ( ! cell_halo) return;
    // For grad v computation, need to collect data in the 1-halo of own cells.
    for (Int j = 0; j < t.vdegree; ++j) {
      const Int cj = t.get_v2c(vi,j);
      if (cj < ncown) continue;
      const Int nv = t.c2ne[cj];
      for (Int k = 0; k < nv; ++k) {
        const Int vk = t.get_c2v(cj,k);
        vtx_mask(li,vk) |= maskval;
        const auto* const stencil = e2v_stencils(vk);
        for (Int l = 0; l < e2v_nedge; ++l) {
          const Int ei = stencil[l];
          if (ei >= neown) {
            assert(ei < nedge);
            edge_mask(li,ei-neown) |= maskval;
          }
        }
        if (is_3d)
          for (Int l = 0; l < t.vdegree; ++l) {
            const Int ci = t.get_v2c(vk,l);
            if (ci >= ncown)
              cell_mask(li,ci-ncown) |= maskval;
          }
      }
    }
  };
  ompfor_ent_lvl(f, l, is_3d, nvown);
}

void Trajectory::Impl
::masks_init_from_cells (const Array2D<const CellLvl>& vtx_in_cell_lvl,
                         const char maskval) {
  const auto& t = mm->t;
  const auto& l = mm->l;
  const auto nvown = mm->nvown, neown = mm->neown, ncown = mm->ncown;
  const auto nedge = t.nedge;
  const auto e2v_nedge = e2v_stencils.m();
  // All of the mask-write lines are benign write race conditions.
  const auto f = [&] (const int vi, const int lie) { // Eulerian inputs
    const auto& cl = vtx_in_cell_lvl(lie,vi);
    const Int ci = cl.cell(), li = cl.level(); // Lagrangian
    assert(ci >= 0 && ci < t.ncell);
    assert( ! is_3d || (li >= 0 && li < l.nmax /* = nilvl-1 */));
    const Int ne = t.c2ne[ci];
    for (Int i = 0; i < ne; ++i) {
      const int vj = t.get_c2v(ci,i);
      assert(vj >= 0 && vj < t.nvtx);
      vtx_mask(li,vj) |= maskval;
      if (is_3d) vtx_mask(li+1,vj) |= maskval;
      const auto* const stencil = e2v_stencils(vj);
      for (Int l = 0; l < e2v_nedge; ++l) {
        const Int ej = stencil[l];
        if (ej >= neown) {
          assert(ej < nedge);
          edge_mask(li,ej-neown) |= maskval;
          if (is_3d) edge_mask(li+1,ej-neown) |= maskval;
        }
      }
      if (is_3d)
        for (Int j = 0; j < t.vdegree; ++j) {
          const Int cj = t.get_v2c(vj,j);
          if (cj >= ncown) {
            cell_mask(li,cj-ncown) |= maskval;
            cell_mask(li+1,cj-ncown) |= maskval;
          }
        }
    }
  };
  ompfor_ent_lvl(f, l, is_3d, nvown);
}

size_t Trajectory::Impl::get_rwrk_size () const {
  const auto& t = mm->t;
  const int dim = is_3d ? 4 : 3;
  return ((size_t(nlvl)*t.nedge) +          // evnmls
          (is_3d ? nlvl*t.ncell : 0) +      // ccz
          (dim*nlvl*t.nvtx) +               // vtxvels
          (method == Method::taylor2 ?      // posits
           nlvl*3*mpas_omp_get_max_threads() :
           method == Method::exact ?
           0 :
           dim*nlvl*mm->nvown));
}

const CellLvl& Trajectory
::get_cell_lvl_containing_vtx (const int vi, const int li) const {
  return impl->vtx_in_cell_lvl(li,vi);
}

Trajectory::OdeFn::OdeFn (const VelocityField& vf_, const int lvl_)
  : vf(vf_), lvl(lvl_), zcomp(vf.is_3d())
{}

bool Trajectory::OdeFn::eval (Real t, const Real* y, Real* f) const {
  const Real zeta = zcomp ? y[2] : 0;
  Real w;
  vf(t, y[0], y[1], zeta, lvl, f[0], f[1], w);
  f[0] /= std::abs(std::cos(y[1]));
  if (zcomp) f[2] = w;
  return true;
}

static Real lininterp (const Real* y, const Real* x, const Real xi) {
  const Real a = (xi - x[0])/(x[1] - x[0]);
  return (1-a)*y[0] + a*y[1];
}

void Trajectory::Impl
::fill_evnmls (const MpasMesh& mm, const Real time, const VelocityField& f,
               const Ar2& ecs, const Acr1& angle_edge, const Ar2& evnmls,
               const int nedge, const int nlvl, Real* const wrk) {
  const bool is_3d = f.is_3d();
  if (is_3d) {
#if 1
    // Edge-normal data live at the level midpoints, so we need to interpolate
    // from there to level interfaces.
    assert(wrk);
    Ar3 awrk(nlvl, 3, mpas_omp_get_max_threads());
    ompfor for (int ei = 0; ei < nedge; ++ei) {
      // Get work arrays.
      const int tid = mpas_omp_get_thread_num();
      Real* const zeta_int = &awrk(0,0,tid);
      Real* const zeta_mid = &awrk(0,1,tid);
      Real* const mevnmls  = &awrk(0,2,tid);
      // Interface zeta.
      const auto& eul_vtx_univalued = mm.zeta_rep->eul_vtx_univalued;
      for (int li = 0; li < nlvl; ++li) {
        Real zeta = 0;
        for (int j = 0; j < 2; ++j) {
          const int vj = mm.t.get_e2v(ei,j);
          zeta += eul_vtx_univalued(li,vj);
        }
        zeta_int[li] = zeta/2;
      }
      // Midpoint zeta and edge velocity normal.
      const Real* const ec = ecs(ei);
      Real lon, lat;
      xyz2lonlat(ec[0], ec[1], ec[2], lon, lat);
      const Real cae = std::cos(angle_edge[ei]), sae = std::sin(angle_edge[ei]);
      for (int li = 0; li < nlvl-1; ++li) {
        zeta_mid[li] = (zeta_int[li] + zeta_int[li+1])/2;
        Real u, v, w;
        f(time, lon, lat, zeta_mid[li], li, u, v, w);
        mevnmls[li] = u*cae + v*sae;
      }
      // Midpoints to interfaces.
      evnmls(0,ei) = lininterp(mevnmls, zeta_mid, zeta_int[0]);
      const int os = nlvl - 3;
      evnmls(nlvl-1,ei) = lininterp(mevnmls + os, zeta_mid + os,
                                    zeta_int[nlvl-1]);
      for (int li = 1; li < nlvl-1; ++li) {
        const int os = li - 1;
        evnmls(li,ei) = lininterp(mevnmls + os, zeta_mid + os, zeta_int[li]);
      }
    }
#else
    // Edge velocity data at level interfaces.
    ompfor for (int ei = 0; ei < nedge; ++ei) {
      const Real* const ec = ecs(ei);
      Real lon, lat;
      xyz2lonlat(ec[0], ec[1], ec[2], lon, lat);
      const Real cae = std::cos(angle_edge[ei]), sae = std::sin(angle_edge[ei]);
      const auto& eul_vtx_univalued = mm.zeta_rep->eul_vtx_univalued;
      for (int li = 0; li < nlvl; ++li) {
        Real zeta = 0;
        for (int j = 0; j < 2; ++j) {
          const int vj = mm.t.get_e2v(ei,j);
          zeta += eul_vtx_univalued(li,vj);
        }
        zeta /= 2;
        Real u, v, w;
        f(time, lon, lat, zeta, li, u, v, w);
        evnmls(li,ei) = u*cae + v*sae;
      }
    }
#endif
  } else {
    ompfor for (int ei = 0; ei < nedge; ++ei) {
      const Real* const ec = ecs(ei);
      Real lon, lat;
      xyz2lonlat(ec[0], ec[1], ec[2], lon, lat);
      const Real cae = std::cos(angle_edge[ei]), sae = std::sin(angle_edge[ei]);
      for (int li = 0; li < nlvl; ++li) {
        Real u, v, w;
        f(time, lon, lat, 0, li, u, v, w);
        evnmls(li,ei) = u*cae + v*sae;
      }
    }    
  }
}

// Fill array of cell center vertical velocities.
void Trajectory::Impl
::fill_ccvverts (const MpasMesh& mm, const Real time, const VelocityField& f,
                 const Ar2& ccs, const Ar2& ccz, const int ncell, const int nlvl) {
  if ( ! f.is_3d()) return;
  ompfor for (int ci = 0; ci < ncell; ++ci) {
    const Real* const cc = ccs(ci);
    Real lon, lat;
    xyz2lonlat(cc[0], cc[1], cc[2], lon, lat);
    const Real* zint = mm.l.get_interface_z(ci);
    for (int li = 0; li < nlvl; ++li) {
      Real u, v, w;
      // Support null zint for testing.
      const Real zeta = zint ? zint[li] : f.ilvl_to_z(nlvl, li);
      f(time, lon, lat, zeta, li, u, v, w);
      ccz(li,ci) = w;
    }
  }
}

void Trajectory::Impl
::get_vtx_v (const Ar2& evnmls, const int vi, const int li, Real vel[3]) const {
  assert(vi >= 0 && vi < e2v_coefs.n());
  const int* stencil = e2v_stencils(vi);
  geo::set(vel, 0);
  const Real* coef = e2v_coefs(vi);
  const int ne = e2v_coefs.m() / 3;
  for (int i = 0; i < ne; ++i) {
    assert(stencil[i] >= 0);
    assert(stencil[i] < evnmls.n());
    geo::axpy(evnmls(li,stencil[i]), &coef[3*i], vel);
  }
}

#undef EXACT_VTX_V

#ifdef EXACT_VTX_V
# pragma message "use get_exact_vtx_v"
static void
get_exact_vtx_v (const MpasMesh& mm, const Real time, const VelocityField& f,
                 const int vi, const int li, Real* const vel) {
  Vec p; fill(mm.v, vi, p, true);
  const Vec axis = {0,0,1};
  Vec east ; geo::cross(axis, p, east ); geo::normalize(east );
  Vec north; geo::cross(p, east, north); geo::normalize(north);
  Real lon, lat, u, v, w;
  xyz2lonlat(p[0], p[1], p[2], lon, lat);
  f(time, lon, lat, 0, li, u, v, w);
  geo::axpbyz(u, east, v, north, vel);
}
#endif

Real Trajectory::Impl
::get_vtx_z (const Ar2& ccz, const int vi, const int li) const {
  assert(vi >= 0 && vi < e2v_coefs.n());
  const auto& t = mm->t;
  Real z = 0;
  int cnt = 0;
  for (int i = 0; i < t.vdegree; ++i) {
    const int ci = t.get_v2c(vi,i);
    if (ci < 0) {
#ifndef NDEBUG
      for (int j = i+1; j < t.vdegree; ++j)
        assert(t.get_v2c(vi,j) < 0);
#endif
      break;
    }
    z += ccz(li,ci);
    ++cnt;
  }
  return z/cnt;
}

template <int dim>
void Trajectory::Impl
::interp_vel_from_vtxs_to1 (const int ci, const int li, const Ar3& vtxvels,
                            const Real* p, Real* p_vel) const {
  const auto& t = mm->t;
  const int nv = t.c2ne[ci];
  const Real* cc = cvsf.ho_coefs(ci);
  const Real* xaxis = cc + 3;
  const Real* yaxis = xaxis + 3;
  const Real* coef_s = yaxis + 3;
  // Fit d component of v by a linear function.
  Real shape[3*dim] = {0};
  for (int j = 0; j < nv; ++j) {
    const Real* const vel = vtxvels(li,t.get_c2v(ci,j));
    for (int i = 0; i < 3; ++i)
      vec<dim>::axpy(coef_s[nv*i + j], vel, &shape[dim*i]);
  }
  // Interpolate to p.
  Vec pmcc; geo::axpbyz(1, p, -1, cc, pmcc);
  const Real x = geo::dot(pmcc, xaxis), y = geo::dot(pmcc, yaxis);
  for (int d = 0; d < dim; ++d)
    p_vel[d] = shape[d] + shape[dim+d]*x + shape[2*dim+d]*y;
}

static void assert_z (const int nlvl, const Real* zint, const Real z) {
#ifndef NDEBUG
  if (!(z >= zint[0     ] - 1e-10*(zint[nlvl-1] - zint[0]))||
      !(z <= zint[nlvl-1] + 1e-10*(zint[nlvl-1] - zint[0])))
    pr(puf(zint[0]) pu(zint[nlvl-1]) pu(z));
  assert(z >= zint[0     ] - 1e-10*(zint[nlvl-1] - zint[0]));
  assert(z <= zint[nlvl-1] + 1e-10*(zint[nlvl-1] - zint[0]));
#endif
}

void Trajectory::Impl
::interp_vel_from_vtxs_to (const CellLvl& cl, const Ar3& vtxvels,
                           const int nilvl, const Real* zint, const Real* p,
                           Real* p_vel) const {
  const auto ci = cl.cell(), li = cl.level();
  if ( ! is_3d) {
    // li is an interface level index.
    interp_vel_from_vtxs_to1<3>(ci, li, vtxvels, p, p_vel);
    return;
  }
  // li is a midpoint level index.
  assert(zint);
  assert(li >= 0 && li < nilvl-1);
  Real z = p[3];
  assert_z(nilvl, zint, z);
  if (z <= zint[0]) z = zint[0];
  else if (z >= zint[nilvl-1]) z = zint[nilvl-1];
  Real v1[4], v2[4];
  interp_vel_from_vtxs_to1<4>(ci, li  , vtxvels, p, v1);
  interp_vel_from_vtxs_to1<4>(ci, li+1, vtxvels, p, v2);
  const Real a = (p[3] - zint[li])/(zint[li+1] - zint[li]);
  vec<4>::axpbyz(1-a, v1, a, v2, p_vel);
}

/* Math for the grad_vel_times_vel* routines:
     p is position (x,y,z) on the unit sphere.
     vel is (xd,yd,zd).
     zeta is vertical coordinate, and zetad is its time derivative.
     dim is 3 for 2D and 4 for 3D.
     v is vel when dim is 3 and (vel,zetad) when dim is 4.
     2D:
       v[d](p) = shape[d] + shape[  dim+d]*xaxis'(p - cc) +
                            shape[2*dim+d]*yaxis'(p - cc)
       v[d]_p = shape[dim+d]*xaxis' + shape[2*dim+d]*yaxis'
       v[d]_p vel = shape[dim+d]*xaxis'vel + shape[2*dim+d]*yaxis'vel
     3D, where v[d,li] refers to the 2D value at interface level li:
       a = (zeta - zint[li])/(zint[li+1] - zint[li])
       a_zeta = 1/(zint[li+1] - zint[li])
       v[d](p,zeta) = (1-a) v[d,li](p) + a v[d,li+1](p)
       v[d]_p = (1-a) v[d,li]_p + a v[d,li+1]_p
       v[d]_zeta = a_zeta (v[d,li+1](p) - v[d,li](p))
       v[d]_(p,zeta) (vel,zetad) = v[d]_p vel + v[d]_zeta zetad
*/

template <int dim>
void Trajectory::Impl
::grad_vel_times_vel1 (const int ci, const int li, const Ar3& vtxvels,
                       const Real* vel, Real* y) {
  const auto& t = mm->t;
  const int nv = t.c2ne[ci];
  const Real* cc = cvsf.ho_coefs(ci);
  const Real* xaxis = cc + 3;
  const Real* yaxis = xaxis + 3;
  const Real* coef_s = yaxis + 3;
  Real shape[2*dim] = {0};
  for (int j = 0; j < nv; ++j) {
    const Real* const vel = vtxvels(li,t.get_c2v(ci,j));
    for (int i = 0; i < 2; ++i)
      vec<dim>::axpy(coef_s[nv*(i+1) + j], vel, &shape[dim*i]);
  }
  const Real xvel = geo::dot(vel, xaxis), yvel = geo::dot(vel, yaxis);
  for (int d = 0; d < dim; ++d)
    y[d] = shape[d]*xvel + shape[dim+d]*yvel;
}

void Trajectory::Impl
::grad_vel_times_vel (const CellLvl& cl, const Ar3& vtxvels,
                      const int nilvl, const Real* zint,
                      const Real* pos, const Real* vel, Real* y) {
  const auto ci = cl.cell(), li = cl.level();
  if ( ! is_3d) {
    // li is an interface level index.
    grad_vel_times_vel1<3>(ci, li, vtxvels, vel, y);
    return;
  }
  // li is a midpoint level index.
  assert(zint);
  assert(li >= 0 && li < nilvl-1);
  Real zeta = pos[3];
  assert_z(nilvl, zint, zeta);
  if (zeta <= zint[0]) zeta = zint[0];
  else if (zeta >= zint[nilvl-1]) zeta = zint[nilvl-1];
  // y = v[d]_p vel.
  const Real a = (zeta - zint[li])/(zint[li+1] - zint[li]);
  Real y1[4], y2[4];
  grad_vel_times_vel1<4>(ci, li  , vtxvels, vel, y1);
  grad_vel_times_vel1<4>(ci, li+1, vtxvels, vel, y2);
  vec<4>::axpbyz(1-a, y1, a, y2, y);
  // y += v[d]_zeta zetad.
  interp_vel_from_vtxs_to1<4>(ci, li  , vtxvels, pos, y1);
  interp_vel_from_vtxs_to1<4>(ci, li+1, vtxvels, pos, y2);
  vec<4>::axpy(-1, y1, y2);
  vec<4>::axpy(vel[3]/(zint[li+1] - zint[li]), y2, y);
}

template <int dim>
void Trajectory::Impl
::grad_vel_times_vel_at_vtx1 (const int vi, const int li, const Ar3& vtxvels,
                              const Real* vel, Real* y) {
  // li is an interface level index.
  const int* vi_v2c = v2c(vi);
  vec<dim>::set(y, 0);
  for (int i = 0; i < vdegree; ++i) {
    Real tmp[dim];
    grad_vel_times_vel1<dim>(vi_v2c[i], li, vtxvels, vel, tmp);
    vec<dim>::axpy(1, tmp, y);
  }
  vec<dim>::scale(1.0/vdegree, y);
}

void Trajectory::Impl
::grad_vel_times_vel_at_vtx (const int vi, const int li, const Ar3& vtxvels,
                             const int nilvl, const Real* vel, Real* y) {
  // li is an interface level index.
  if ( ! is_3d) {
    grad_vel_times_vel_at_vtx1<3>(vi, li, vtxvels, vel, y);
    return;
  }
  if (li == 0 || li == nilvl-1) {
    grad_vel_times_vel_at_vtx1<4>(vi, li, vtxvels, vel, y);
    return;
  }
  //todo-perf Version that doesn't compute unneeded components 0:2 twice.
  Real y1[4], y2[4];
  grad_vel_times_vel_at_vtx1<4>(vi, li-1, vtxvels, vel, y1);
  grad_vel_times_vel_at_vtx1<4>(vi, li  , vtxvels, vel, y2);
  vec<4>::axpbyz(0.5, y1, 0.5, y2, y);
}

// We can use the Voronoi property. p is in the cell whose cell center is
// nearest. Can use Euclidean distance, too.
//todo-pebfb Settle ties using cell GID.
inline int Trajectory::Impl::get_containing_cell_given_relations (
  const int* const relation, const int nrelation, const Real* const p)
{
  Real min_dist2 = 1e30;
  int min_ci = -1;
  for (int i = 0; i < nrelation; ++i) {
    const int ci = relation[i];
    assert(ci >= 0);
    const Real* cc = cvsf.ho_coefs(ci);
    Vec pmcc; geo::axpbyz(1, p, -1, cc, pmcc);
    const Real d2 = geo::norm2(pmcc);
    if (d2 < min_dist2) {
      min_dist2 = d2;
      min_ci = ci;
    }
  }
  assert(min_ci >= 0);
  return min_ci;
}

#undef TRAJ_OLD_SEARCH

int Trajectory::Impl
::get_containing_cell_given_vtx (const int vi, const Real* p) {
#ifdef TRAJ_OLD_SEARCH
  int relations[(3*max_6halo_slots)/2];
  const int nnbr = v2ccnt[vi];
  assert(nnbr > 0);
  const int* vi_v2c = v2c(vi);
  for (int i = 0; i < nnbr; ++i) relations[i] = vi_v2c[i];
  return get_containing_cell_given_relations(relations, nnbr, p);
#else
  return find_cell_given_starting_cell(mm->t.get_v2c(vi,0), p);
#endif
}

int Trajectory::Impl
::get_containing_cell_given_cell (const int ci, const Real* p) {
#ifdef TRAJ_OLD_SEARCH
  int relations[max_6halo_slots];
  const int nnbr = c2ccnt[ci];
  relations[0] = ci;
  if (nnbr <= 0) return -1;
  const int* ci_c2c = c2c(ci);
  for (int i = 0; i < nnbr; ++i) relations[i+1] = ci_c2c[i];
  return get_containing_cell_given_relations(relations, nnbr+1, p);
#else
  return find_cell_given_starting_cell(ci, p);
#endif
}

// Return the level (midpoint index) containing z.
int Trajectory::Impl::get_containing_level_given_z (
  const int nilvl,        // number level interfaces
  const Real* const zint, // z interface values; zint[k+1] >= zint[k]
  const Real z,           // find midpoint lvl idx for this value
  const int lvl_hint)     // optionally start here
{
  assert(nilvl >= 2);
  if (lvl_hint >= 0 && lvl_hint < nilvl-1 &&
      z >= zint[lvl_hint] && z <= zint[lvl_hint+1])
    return lvl_hint;
  assert_z(nilvl, zint, z);
  if (z <= zint[0]) return 0;
  if (z >= zint[nilvl-1]) return nilvl-2;
  const int idx = std::lower_bound(zint, zint+nilvl, z) - zint - 1;
  assert(idx >= 0 && idx < nilvl-1);
  return idx;
}

CellLvl Trajectory::Impl
::get_containing_cell_lvl_given_vtx_lvl (const Levels& lmgr, const int vi,
                                         const int li, const Real* p) {
  const int ci = get_containing_cell_given_vtx(vi, p);
  if (is_3d) {
    const int li_new = get_containing_level_given_z(
      lmgr.nmax+1, lmgr.get_interface_z(ci), p[3], li);
    return CellLvl(ci, li_new);
  } else {
    return CellLvl(ci, li);  
  }
}

CellLvl Trajectory::Impl
::get_containing_cell_lvl_given_cell_lvl (const Levels& lmgr, const CellLvl& cl,
                                          const Real* p) {
  const int ci = get_containing_cell_given_cell(cl.cell(), p);
  if (is_3d) {
    const int li = get_containing_level_given_z(
      lmgr.nmax+1, lmgr.get_interface_z(ci), p[3], cl.level());
    return CellLvl(ci, li);
  } else {
    return CellLvl(ci, cl.level());
  }
}

static void calc_step (const bool is_3d, const Real dt, const Real* vel, Real* p) {
#if 1
  geo::axpy(dt, vel, p);
  geo::normalize(p);
#else
  // Correct for sphere.
  Vec p0; geo::copy(p0, p);
  geo::axpy(dt, vel, p);
  Vec diff; geo::axpbyz(1, p, -1, p0, diff);
  const Real len = geo::norm(diff);
  geo::normalize(p);
  for (int it = 0; it < 1; ++it) {
    const Real scale = len / geo::calc_arc_length(p, p0);
    geo::scale(scale, diff);
    geo::axpbyz(1, p0, 1, diff, p);
    geo::normalize(p);
  }
#endif
  if (is_3d) p[3] += dt*vel[3];
}

struct ErrorReport {
  bool error;
  Int vi, li, si, it;
  Vec p;
  ErrorReport () : error(false) {}
};

static void set_error (int vi, int li, int si, int it, const Vec p, ErrorReport& er) {
  er.error = true;
  er.vi = vi; er.li = li; er.si = si; er.it = it;
  for (Int d = 0; d < 3; ++d) er.p[d] = p[d];
}

static void report_error_if_needed (const ErrorReport& er) {
  require_msg( ! er.error,
               "Trajectory: Lagrangian vertex, level " << er.vi << " " << er.li
               << " out of bounds at si,it " << er.si << "," << er.it
               << " with p (" << er.p[0] << "," << er.p[1] << "," << er.p[2]
               <<  ")");
}

//todo-perf template on dim and then make vec<dim> calls instead of
//          geo::*n(dim,...) ones.
void Trajectory::Impl
::advect_vtxs (const Acr2& p, const Real ts, const Real dt, const Ar1& wrk,
               const Ar3& p_advect) {
  // When is_3d, pos and vel have components (x,y,z,zeta), where (x,y,z) is
  // position on the unit circle and zeta is vertical coordinate.
  const auto& t = mm->t;
  Ar2 evnmls, ccz;
  Ar3 vtxvels, posits;
  auto* pwrk = wrk.data();
  evnmls.reset(nlvl, t.nedge, pwrk); pwrk += evnmls.size();
  if (is_3d) { ccz.reset(nlvl, t.ncell, pwrk); pwrk += ccz.size(); }
  vtxvels.reset(is_3d ? 4 : 3, nlvl, t.nvtx, pwrk); pwrk += vtxvels.size();
  const auto& lmgr = mm->l;
  const auto& dist = he->cdhalo();
  Ar3 pevnmls_f(1, nlvl, mm->t.nedge, evnmls.data());
  Ar3 ccv_f(1, nlvl, mm->t.ncell, ccz.data());
  const int vflds_ns[] = {1};
  assert(p_advect.m() == 4 || ! is_3d);
  assert(p_advect.n() >= nlvl);
  const int dim = is_3d ? 4 : 3;
  const auto& ccs = mm->g.my_cellCentroid;
  ErrorReport er;
  switch (method) {
  case Method::exact: {
    for (Int li = 0; li < nlvl; ++li) {
      const OdeFn f_ode(*vf, li);
      const auto& eul_vtx_univalued = mm->zeta_rep->eul_vtx_univalued;
      ompparfor for (int vi = 0; vi < mm->nvown; ++vi) {
        timeint::Options opts;
        opts.set_abs_tol(std::numeric_limits<Real>::epsilon());
        opts.set_rel_tol(method_exact_rel_tol);
        opts.set_initial_step(initial_step[vi]);
        Real pos[3], wrk[27];
        xyz2lonlat(p(0,vi), p(1,vi), p(2,vi), pos[0], pos[1]);
        pos[2] = is_3d ? eul_vtx_univalued(li,vi) : 0;
        timeint::Info info;
        timeint::ark45(opts, f_ode, pos, is_3d ? 3 : 2, ts + dt, ts, wrk, &info);
        const auto pa = p_advect(li,vi);
        lonlat2xyz(pos[0], pos[1], pa[0], pa[1], pa[2]);
        if (is_3d) pa[3] = pos[2];
        if (li == nlvl-1) initial_step[vi] = info.good_initial_step;
      }
    }
  } break;
  case Method::order1:
  case Method::order2: {
    // Simple iterative solver for midpoint velocity that is consistent with
    // trajectory. One iteration gives OOA 1, and >= 2 give OOA 2.
    const int nit = method == Method::order1 ? 1 : 2;
    if (posits.empty()) {
      posits.reset(p_advect.m(), nlvl, mm->nvown, pwrk);
      pwrk += posits.size();
      assert( ! is_3d || (posits.size() >= size_t(nlvl)*3*mpas_omp_get_max_threads()));
    }
    const Real dtsub = dt/nsubstep;
    ompparfor for (int vi = 0; vi < mm->nvown; ++vi) {
      for (Int li = 0; li < nlvl; ++li)
        geo::copy(p_advect(li,vi), p(vi));
      if (is_3d) {
        const auto zint = mm->zeta_rep->eul_vtx_univalued(vi);
        for (Int li = 0; li < nlvl; ++li)
          p_advect(3,li,vi) = zint[li];
      }
    }
    for (int si = 0; si < nsubstep; ++si) {
      omppar {
        const auto time = ts + dt - (si+0.5)*dtsub;
        fill_evnmls(*mm, time, *vf, ecs, angle_edge, evnmls, mm->neown, nlvl,
                    is_3d ? posits.data() : nullptr);
        fill_ccvverts(*mm, time, *vf, ccs, ccz, mm->ncown, nlvl);
      }
      for (int it = 0; it < nit; ++it) {
        const char maskval =
          (char) (it == 0 && si == 0 ? MaskVals::from_vtx :
                  it == 0 ? MaskVals::from_cell : MaskVals::from_cell_it1);
        if (si > 0 || it > 0) {
          omppar {
            // In its > 0, don't recommunicate or recompute values that we
            // already have from it 0.
            if (it == 0) {
              mask_and(edge_mask, (char) MaskVals::from_vtx);
              mask_and(vtx_mask,  (char) MaskVals::from_vtx);
              mask_and(cell_mask, (char) MaskVals::from_vtx);
            }
            masks_init_from_cells(vtx_in_cell_lvl, maskval);
          }
        }
        if (is_3d)
          he->exchange_edge_and_cell_masked(edge_mask, cell_mask, maskval, 1,
                                            vflds_ns, nlvl, &pevnmls_f, &ccv_f);
        else
          he->exchange_edge_masked(edge_mask, maskval, 1, vflds_ns, nlvl,
                                   &pevnmls_f);
        omppar {
          const auto f = [&] (const int vi, const int li) {
            if ( ! (vtx_mask(li,vi) & maskval)) return;
            get_vtx_v(evnmls, vi, li, vtxvels(li,vi));
            if (is_3d) vtxvels(3,li,vi) = get_vtx_z(ccz, vi, li);
          };
          ompfor_ent_lvl(f, mm->l, is_3d, t.nvtx);
          const auto g = [&] (const int vi, const int li) {
            Real* const p_a = p_advect(li,vi);
            Real vel[4] = {0};
            Real* const pos = posits(li,vi);
            const auto cl_cur = vtx_in_cell_lvl(li,vi);
            if (it == 0) geo::copyn(dim, pos, p_a);
            if (si == 0 && it == 0)
              geo::copyn(dim, vel, vtxvels(li,vi));
            else {
              const Real* zint = mm->l.get_interface_z(cl_cur.cell());
              interp_vel_from_vtxs_to(cl_cur, vtxvels, nlvl, zint, pos, vel);
            }
            Real* p;
            bool use_vtx;
            if (it+1 < nit) {
              geo::axpbyzn(dim, -0.5*dtsub, vel, 1, p_a, pos);
              p = pos;
              use_vtx = si == 0 && it == 0;
            } else {
              geo::axpyn(dim, -dtsub, vel, p_a);
              p = p_a;
              use_vtx = si == 0 && nit == 1;
            }
            geo::normalize(p); // 3-vec
            const auto cl_new = 
              use_vtx ?
              get_containing_cell_lvl_given_vtx_lvl(lmgr, vi, li, p) :
              get_containing_cell_lvl_given_cell_lvl(lmgr, cl_cur, p);
            vtx_in_cell_lvl(li,vi) = cl_new;
            const auto ci = cl_new.cell();
            assert(ci >= 0 && ci < mm->t.ncell);
            if (ci < 0 || dist(ci) == 1)
              set_error(vi, li, si, it, it+1 < nit ? pos : p_a, er);
          };
          ompfor_ent_lvl(g, mm->l, is_3d, mm->nvown);
        }
        if (er.error) break;
      }
      if (er.error) break;
    }
  } break;
  case Method::taylor2: {
    // One iteration using grad v to boost OOA to 2.
    //     x0 = x1 - dt v(xm,tm) + O(dt^3)
    //     v(xm,tm) = v(x1,tm) + grad v(x1,tm) (xm - x1) + O(dt^2)
    //     xm - x1 = -dt/2 v(x1,tm) + O(dt)
    //     x0 = x1 - dt [v(x1,tm) - dt/2 grad v(x1,tm) v(x1,t1)] + O(dt^3)
    if (is_3d && posits.empty()) {
      posits.reset(nlvl, 3, mpas_omp_get_max_threads(), pwrk);
      pwrk += posits.size();
    }
    const Real dtsub = dt/nsubstep;
    for (int si = 0; si < nsubstep; ++si) {
      omppar {
        const auto time = ts + dt - (si+0.5)*dtsub;
        fill_evnmls(*mm, time, *vf, ecs, angle_edge, evnmls, mm->neown, nlvl,
                    is_3d ? posits.data() : nullptr);
        fill_ccvverts(*mm, time, *vf, ccs, ccz, mm->ncown, nlvl);
        if (si > 0) {
          // Zero out the from_cell bit in eash vertex mask.
          mask_and(edge_mask, (char) MaskVals::from_vtx);
          mask_and(vtx_mask,  (char) MaskVals::from_vtx);
          mask_and(cell_mask, (char) MaskVals::from_vtx);
          // Set the from_cell bit for each vertex whose velocity is needed.
          masks_init_from_cells(vtx_in_cell_lvl, (char) MaskVals::from_cell);
        }
      }
      const char maskval = (char) (si == 0 ?
                                   MaskVals::from_vtx :  // set once at init
                                   MaskVals::from_cell); // set above
      if (is_3d)
        he->exchange_edge_and_cell_masked(edge_mask, cell_mask, maskval, 1,
                                          vflds_ns, nlvl, &pevnmls_f, &ccv_f);
      else
        he->exchange_edge_masked(edge_mask, maskval, 1, vflds_ns, nlvl,
                                 &pevnmls_f);
      omppar {
        const auto f = [&] (const int vi, const int li) {
          if ( ! (vtx_mask(li,vi) & maskval)) return;
          // Needed if not using masked comm to omit the shared perimeter.
          //   if (v2ccnt[vi] == 0) return;
          get_vtx_v(evnmls, vi, li, vtxvels(li,vi));
          if (is_3d) vtxvels(3,li,vi) = get_vtx_z(ccz, vi, li);
        };
        ompfor_ent_lvl(f, mm->l, is_3d, t.nvtx);
        const auto g = [&] (const int vi, const int li) {
          Real* p_a = p_advect(li,vi);
          if (si == 0) {
            geo::copy(p_a, p(vi));
            if (is_3d) p_a[3] = mm->zeta_rep->eul_vtx_univalued(li,vi);
          }
          const auto cl_cur = vtx_in_cell_lvl(li,vi);
          Real vel[4], cterm[4];
          if (si == 0) {
            geo::copyn(dim, vel, vtxvels(li,vi));
            grad_vel_times_vel_at_vtx(vi, li, vtxvels, nlvl, vel, cterm);
          } else {
            const Real* zint = mm->l.get_interface_z(cl_cur.cell());
            interp_vel_from_vtxs_to(cl_cur, vtxvels, nlvl, zint, p_a, vel);
            grad_vel_times_vel(cl_cur, vtxvels, nlvl, zint, p_a, vel, cterm);
          }
          geo::axpyn(dim, -0.5*dtsub, cterm, vel);
          calc_step(is_3d, -dtsub, vel, p_a);
          const auto cl_new =
            si == 0 ?
            get_containing_cell_lvl_given_vtx_lvl(lmgr, vi, li, p_a) :
            get_containing_cell_lvl_given_cell_lvl(lmgr, cl_cur, p_a);
          vtx_in_cell_lvl(li,vi) = cl_new;
          const auto ci = cl_new.cell();
          if (ci < 0 || dist(ci) == 1) set_error(vi, li, si, 0, p_a, er);
        };
        ompfor_ent_lvl(g, mm->l, is_3d, mm->nvown);
      }
      if (er.error) break;
    }
  } break;
  default: assert(0);
  }
  report_error_if_needed(er);
  assert(pwrk >= wrk.data() && size_t(pwrk - wrk.data()) <= get_rwrk_size());
  assert(size_t(pwrk - wrk.data()) <= wrk.size());
}

size_t Trajectory::get_rwrk_size () const { return impl->get_rwrk_size(); }

bool Trajectory::has_in_cell_data () const {
  return ! Method::is_exact(impl->method);
}

void Trajectory
::advect_vtxs (const Acr2& p, const Real ts, const Real dt, const Ar1& rwrk,
               const Ar3& p_advect) {
  impl->advect_vtxs(p, ts, dt, rwrk, p_advect);
}

void adjust_nsubstep_nhalo_if_needed (
  const Comm& he, const Trajectory::Method::Enum method, const int delta,
  int& nsubstep, int& nhalo)
{
  if (( ! he.parallel() && he.szhalo() == 0) ||
      Trajectory::Method::is_exact(method))
    return;
  // Assume the comm halo ~ the max trajectory length in a step.
  const int nhalo_needed = (nsubstep == 0 ?
                            max_nhalo :
                            (he.szhalo()-delta + nsubstep - 1)/nsubstep);
  if (nhalo < nhalo_needed)
    nhalo = std::min(max_nhalo, nhalo_needed);
  if (nhalo_needed > max_nhalo)
    nsubstep = std::max(nsubstep,
                        (he.szhalo()-delta + nhalo - 1)/nhalo);
}
