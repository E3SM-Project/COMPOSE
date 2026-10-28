#include "mpas_trajectory_impl.hpp"

struct Trajectory::Impl::TestZVelocityField : public VelocityField {
  bool is_3d () const override { return true; }

  void operator() (const Real, const Real lon, const Real lat, const Real,
                   const Int, Real& u, Real& v, Real& w) const override {
    u = v = 0;
    Real x, y, z;
    lonlat2xyz(lon, lat, x, y, z);
    w = 0.5*(1 + std::sin(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z));
  }

  static bool test_z_calcs (Trajectory::Impl& p) {
    bool ok = true;
    const int nvtx = p.mm->t.nvtx, ncell = p.mm->t.ncell, nlvl = 2;
    Ar2 ccz(nlvl, ncell);
    Ar3 vtxvels(4, nlvl, nvtx);
    TestZVelocityField f;
    omppar {
      fill_ccvverts(*p.mm, 0, f, p.mm->g.my_cellCentroid, ccz, ncell, nlvl);
      ompfor for (int vi = 0; vi < nvtx; ++vi)
        for (int li = 0; li < nlvl; ++li)
          vtxvels(3,li,vi) = p.get_vtx_z(ccz, vi, li);
    }
    const int li = 0;
    Real z_num = 0, z_den = 0;
    const Real zint[] = {0, 1}, z = 0.5;
#   pragma omp parallel reduction (+: z_num, z_den)
    for (int ci = 0; ci < ncell; ++ci) {
      const auto area = p.mm->g.areaCell[ci];
      const auto cc = p.mm->g.my_cellCentroid(ci);
      Real pos[4], vel[4];
      geo::copy(pos, cc);
      pos[3] = z;
      p.interp_vel_from_vtxs_to(CellLvl(ci, li), vtxvels, nlvl, zint, pos, vel);
      Real lon, lat, ut, vt, wt;
      xyz2lonlat(cc[0], cc[1], cc[2], lon, lat);
      f(0, lon, lat, 0, 0, ut, vt, wt);
      z_num += area*square(wt - vel[3]);
      z_den += area*square(wt);
    }
    const Real e = std::sqrt(z_num/z_den);
    printf("re l2 z %1.3e\n", e);
    if (e > 1e-2) ok = false;
    return ok;
  }
};

bool Trajectory::Impl
::unittest (const Acr2& vtxs, Real* data) {
  bool ok = true;
  const int nvtx = mm->t.nvtx, ncell = mm->t.ncell;

  // 2D tests.
  bool is_3d_orig = is_3d;
  is_3d = false;

  { // Test augment_stencil.
    int nerr = 0;
    std::vector<int> s1(128), s2(128);
    const int nhalo = 6;
    const int ci = mm->t.ncell / 2;
    int cnt1 = 0, cnt2;
    for (int halo = 0; halo < nhalo; ++halo) {
      cnt1 = augment_stencil(*mm, ci, cnt1, s1.data(), 1);
      printf("halo %d cnt %d\n", halo+1, cnt1+1);
    }
    cnt2 = augment_stencil(*mm, ci, 0, s2.data(), nhalo);
    if (cnt1 != cnt2) {
      printf("FAIL augment_stencil cnt1 == cnt2: %d vs %d\n", cnt1, cnt2);
      ++nerr;
    } else {
      for (int i = 0; i < cnt1; ++i) if (s1[i] != s2[i]) ++nerr;
      if (nerr) printf("FAIL augment_stencil cmp: cnt %d nerr %d\n", cnt1, nerr);
      printf("cnt %d\n",cnt1);
    }
    if (nerr) ok = false;
  }

  { // Test get_containing_cell*.
    const int nrand = 3;
    int nerr = 0;
    const auto test_cell = [&] (const bool iscell, const int idx, const int ci) {
      int my_nerr = 0;
      Real a[max_1halo_slots];
      const int ne = mm->t.c2ne[ci];
      // Generate a random point inside this cell.
      for (int i = 0; i < ne; ++i) a[i] = urand();
      Real sum = 0; for (int i = 0; i < ne; ++i) sum += a[i];
      for (int i = 0; i < ne; ++i) a[i] /= sum;
      Vec p = {0};
      for (int i = 0; i < ne; ++i)
        geo::axpy(a[i], vtxs(mm->t.get_c2v(ci,i)), p);
      geo::normalize(p);
      const int ci_fnd = (iscell ?
                          get_containing_cell_given_cell(idx, p) :
                          get_containing_cell_given_vtx(idx, p));
      if (ci_fnd != ci) ++my_nerr;
      return my_nerr;
    };
#   pragma omp parallel reduction (+: nerr)
    for (int vi = 0; vi < nvtx; ++vi)
      for (int i = 0; i < vdegree; ++i)
        for (int k = 0; k < nrand; ++k)
          nerr += test_cell(false, vi, mm->t.get_v2c(vi,i));
#   pragma omp parallel reduction (+: nerr)
    for (int ci = 0; ci < ncell; ++ci) {
      const int ne = mm->t.c2ne[ci];
      const int* ci_c2c = c2c(ci);
      for (int k = 0; k < nrand; ++k)
        nerr += test_cell(true, ci, ci);
      for (int i = 0; i < ne; ++i)
        for (int k = 0; k < nrand; ++k)
          nerr += test_cell(true, ci, ci_c2c[i]);
    }
    if (nerr) printf("get_containing_cell FAILed %d times\n", nerr);
    if (nerr) ok = false;
  }

  { // Test get_vtx_v and interp_vel_from_vtxs_to.
    const Real time = 3.5*86400;
    Ar2 evnmls(nlvl, mm->t.nedge);
    Ar3 vtxvels(3, nlvl, nvtx);
    omppar {
      fill_evnmls(*mm, time, *vf, ecs, angle_edge, evnmls, mm->t.nedge, nlvl);
      ompfor for (int vi = 0; vi < nvtx; ++vi)
        for (int li = 0; li < nlvl; ++li) {
          Real* const vel = vtxvels(li,vi);
          get_vtx_v(evnmls, vi, li, vel);
        }
    }
    Real u_num = 0, u_den = 0, v_num = 0, v_den = 0;
#   pragma omp parallel reduction (+: u_num, u_den, v_num, v_den)
    for (int ci = 0; ci < ncell; ++ci) {
      const auto area = mm->g.areaCell[ci];
      for (int li = 0; li < nlvl; ++li) {
        const Real* cc = cvsf.ho_coefs(ci);
        Real vel[3];
        interp_vel_from_vtxs_to(CellLvl(ci, li), vtxvels, 0, nullptr, cc, vel);
        Real lon, lat, u, v, ut, vt, wt;
        xyz2lonlat(cc[0], cc[1], cc[2], lon, lat);
        Vec axis = {0,0,1}, east, nrth;
        if (std::abs(lat) < (1 - 1e-8)*0.5*M_PI) {
          geo::cross(axis, cc, east); geo::normalize(east);
          geo::cross(cc, east, nrth); geo::normalize(nrth);
          u = geo::dot(vel, east);
          v = geo::dot(vel, nrth);
        } else {
          u = v = 0;
        }
        (*vf)(time, lon, lat, 0, li, ut, vt, wt);
        u_num += area*square(ut - u); v_num += area*square(vt - v);
        u_den += area*square(ut); v_den += area*square(vt);
        if (data) {
          data[0*ncell + ci] = u;  data[1*ncell + ci] = v;
          data[2*ncell + ci] = ut; data[3*ncell + ci] = vt;
        }
      }
    }
    const Real e1 = std::sqrt(u_num/u_den), e2 = std::sqrt(v_num/v_den);
    printf("re l2 u %1.3e v %1.3e\n", e1, e2);
    if (std::max(e1, e2) > 1e-2) ok = false;
  }

  // 3D tests.
  is_3d = true;
  
  // Test get_vtx_z and interp_z_from_vtxs_to.
  if ( ! TestZVelocityField::test_z_calcs(*this)) ok = false;

  { // Test get_containing_level_given_z.
    const Real zint[] = {0.5, 1, 10, 11.1};
    const int nilvl = sizeof(zint)/sizeof(*zint);
    const auto test = [&] (Real z, int lvl, int hint) {
      const int mlvl = get_containing_level_given_z(nilvl, zint, z, hint);
      if (mlvl != lvl) {
        ok = false;
        pr("get_containing_level_given_z test FAIL:"
           pu(nilvl) pu(z) pu(hint) pu(lvl) pu(mlvl));
      }
    };
    test(zint[0], 0, -1);
    test(zint[nilvl-1], nilvl-2, -1);
    test(0.7, 0, -1);
    test(0.7, 0, 12);
    test(10.1, 2, 12);
  }

  is_3d = is_3d_orig;
  return ok;
}

static Real reldif (const Acr2& a, const Array3D<const Real>& b, const int li) {
  Real num = 0, den = 0;
  const int n = nslices(a);
  for (int j = 0; j < n; ++j) {
    Vec d; geo::axpbyz(1, a(j), -1, b(li,j), d);
    num += geo::norm2(d);
    den += geo::norm2(a(j));
  }
  return std::sqrt(num/den);
}

bool Trajectory::Impl
::convtest (const MpasMesh::ConstPtr& m, const Comm::Ptr& he, Real* data) {
  // This test is meant to be used with a sequence of grids, but we'll do some
  // basic checks, too.
  bool ok = true;
  
  const auto vf = vf_from_short('n');
  const Int nlvl = m->l.nmax + (vf->is_3d() ? 1 : 0);

  Array2D<Real> p0(3, m->t.nvtx);
  ompparfor for (int i = 0; i < m->t.nvtx; ++i) {
    p0(0,i) = m->v.x[i]; p0(1,i) = m->v.y[i]; p0(2,i) = m->v.z[i];
    geo::normalize(p0(i));
  }

  Real R, D;
  get_avg_radius_and_cell_spacing(*m, *he, R, D);
  const Real T = vf->cycle_time(), max_vel = vf->max_vel(), cfl = (D/R)/max_vel;
  const int base_nsubstep = (int) (T/cfl + 1);

  Ar3 pf(3, nlvl, m->t.nvtx);
  {
    Trajectory te(m, he, vf, Trajectory::Method::exact, -1);
    Ar1 rwrk(te.get_rwrk_size());
    te.advect_vtxs(p0, T, T, rwrk, pf);
  }
  for (int li = 0; li < nlvl; ++li) {
    const auto e = reldif(p0, pf, li);
    printf("exact %1.3e\n", e);
    if (e > 1e-6) ok = false;
  }
  printf("\n");
  // These next methods depend on operations in space, so convergence cannot be
  // assessed from a single grid. Instead, take numbers from a sequence of calls
  // to convtest on a sequence of grids.
  struct {
    Trajectory::Method::Enum method;
    const char* name;
  } methods[] = {{Trajectory::Method::order1, "o1"},
                 {Trajectory::Method::order2, "o2"},
                 {Trajectory::Method::taylor2, "t2"}};
  for (size_t mi = 0; mi < sizeof(methods)/sizeof(*methods); ++mi) {
    const auto& method = methods[mi];
    for (int i = 0; i < 4; ++i) {
      const int nsubstep = base_nsubstep*(1 << i);
      Trajectory traj(m, he, vf, method.method, nsubstep);
      Ar1 rwrk(traj.get_rwrk_size());
      traj.advect_vtxs(p0, T, T, rwrk, pf);
      for (int li = 0; li < nlvl; ++li) {
        const auto e = reldif(p0, pf, li);
        printf("%s %d %1.3e\n", method.name, i, e);
        if (i == 3 && e > 0.1) ok = false;
      }
    }
    printf("\n");
  }
  if ( ! ok) printf("Trajectory::convtest FAILed\n");
  return ok;
}

bool Trajectory
::unittest (const Acr2& vtxs, Real* data) {
  return impl->unittest(vtxs, data);
}

bool Trajectory
::convtest (const MpasMesh::ConstPtr& m, const Comm::Ptr& he, Real* data) {
  return Impl::convtest(m, he, data);
}

bool unittest_gallery (const MpasMesh::Ptr& mm, const siqp::sh::Mesh& m,
                       const Comm::Ptr& he) {
  bool ok = true;
  const Int nvtx = mm->nvown;
  for (const char vf_short : {'n', '3', 'd', 'r', 'a'}) {
    const auto vf = vf_from_short(vf_short);
    set_interface_z_all_cols(*vf, *mm);
    init_zeta_rep(mm, m, he);
    const Int nlvl = mm->l.nmax + (vf->is_3d() ? 1 : 0);
    Ar3 p_advect(4, nlvl, nvtx);
    const auto H = std::abs(vf->ilvl_to_z(2, 0) - vf->ilvl_to_z(2, 1));
    const bool zcomp = H > 0;
    assert(vf->is_3d() == zcomp);
    const auto traj = std::make_shared<Trajectory>(mm, he, vf,
                                                   Trajectory::Method::exact, 1);
    Ar1 rwrk(traj->get_rwrk_size());
    const Real reltol = 1e6*std::numeric_limits<Real>::epsilon();
    traj->set_method_exact_rel_tol(reltol);
    traj->advect_vtxs(m.p, 0, vf->cycle_time(), rwrk, p_advect);
    // Check that after one cycle, each point returns to its starting point to
    // w/in the tolerance of the ODE integrator.
    Real max_err_sphere = 0, max_err_z = 0;
    for (int vi = 0; vi < nvtx; ++vi)
      for (int li = 0; li < nlvl; ++li) {
        Real err = 0;
        for (int d = 0; d < 3; ++d)
          err += square(p_advect(d,li,vi) - m.p(d,vi));
        max_err_sphere = std::max(max_err_sphere, err); // w.r.t. unit sphere
        max_err_z = std::max(max_err_z, std::abs(p_advect(3,li,vi) -
                                                 vf->ilvl_to_z(nlvl, li)));
      }
    max_err_sphere = std::sqrt(max_err_sphere);
    if (zcomp) max_err_z /= H;
    const Real l[2] = {max_err_sphere, max_err_z};
    Real g[2];
    he->all_reduce(Comm::MAX, 2, l, g);
    if (g[0] > (zcomp? 50 : 10)*reltol) {
      ok = false;
      if (he->amroot())
        printf("unittest_gallery %c sphere relerr %1.3e\n", vf_short, g[0]);
    }
    if (zcomp && g[1] > 10*reltol) {
      ok = false;
      if (he->amroot())
        printf("unittest_gallery %c z relerr %1.3e\n", vf_short, g[1]);
    }
  }
  return ok;
}
