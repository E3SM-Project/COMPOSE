#include "mpas_check.hpp"
#include "mpas_geometry.hpp"
#include "mpas_advect_ir_test.hpp"
#include "mpas_siqp_test.hpp"
#include "mpas_gallery.hpp"
#include "mpas_trajectory.hpp"

#include <set>
#include <vector>
#include <algorithm>

extern "C" void unit_test () {
  if (timeint::test::test_ark(false) != 0) {
    printf("FAIL test_ark\n");
    timeint::test::test_ark(true);
  }
}

static MpasMesh::Ptr s_mm;
static siqp::sh::Mesh::Ptr s_m;

extern "C" void set_mesh (
  int ncell, int nedge, int nvtx, int max_c2e, int nvdegree,
  const PyInt* c2ne, const PyInt* c2c, const PyInt* c2e, const PyInt* c2v,
  const PyInt* e2c, const PyInt* e2v,
  const PyInt* v2c, const PyInt* v2e,
  const Real* xc, const Real* yc, const Real* zc,
  const Real* xe, const Real* ye, const Real* ze,
  const Real* xv, const Real* yv, const Real* zv,
  const Real* dvEdge, const Real* dcEdge, const Real* angleEdge, const Real* areaCell,
  // For edge velocity reconstruction
  const PyInt* e2e, const PyInt* e2ne, const Real* e2w, int max_e2w)
{
  s_mm = std::make_shared<MpasMesh>(
    ncell, nedge, nvtx, max_c2e, nvdegree,
    c2ne, c2c, c2e, c2v, e2c, e2v, v2c, v2e,
    xc, yc, zc, xe, ye, ze, xv, yv, zv,
    dvEdge, dcEdge, angleEdge, areaCell,
    e2e, e2ne, e2w, max_e2w);
  init_extras(*s_mm);
  s_m = create_siqp_mesh(*s_mm);
  init_more_extras(*s_mm, *s_m);
}

extern "C" void test_siqp () {
  assert(s_mm);
  const auto he = Comm::create();
  if ( ! test_siqp(he, s_mm)) printf("FAIL test_siqp\n");
}

extern "C" void check (const int out_opt, Real* out) {
  assert(s_mm);
  check(*s_mm);
  if (out_opt == 42) return;
  switch (out_opt) {
  case 0: get_angleEdge_error(s_mm->t, s_mm->c, s_mm->e, s_mm->v, s_mm->g, out); break;
  default: break;
  }
}

extern "C" void get_oriented_cells (int* elems) {
  assert(s_mm);
  for (int ci = 0; ci < s_mm->t.ncell; ++ci)
    get_oriented_cell(s_mm->t, s_mm->v, ci, elems + s_mm->t.max_vtx * ci);
}

static void wrap_lonlat (Real* xs, Real* ys, const int ne) {
  Real min_lon, max_lon, min_lat, max_lat;
  const auto calc_ext = [&] () {
    min_lon = 1e30; max_lon = -1e30; min_lat = 1e30; max_lat = -1e30;
    for (int i = 0; i < ne; ++i) {
      min_lon = std::min(min_lon, xs[i]);
      max_lon = std::max(max_lon, xs[i]);
      min_lat = std::min(min_lat, ys[i]);
      max_lat = std::max(max_lat, ys[i]);
    }      
  };
  calc_ext();
  if (max_lat < -0.95*M_PI/2 || min_lat > 0.95*M_PI/2) return;
  if (max_lat > 0 && min_lat < 0 && max_lat - min_lat > M_PI/2)
    for (int i = 0; i < ne; ++i)
      if (ys[i] < 0) ys[i] += M_PI;
  if (max_lon > 0 && min_lon < 0 && max_lon - min_lon > M_PI) {
    if (std::abs(max_lon) > std::abs(min_lon)) {
      for (int i = 0; i < ne; ++i)
        if (xs[i] > M_PI/2) xs[i] -= 2*M_PI;
    } else {
      for (int i = 0; i < ne; ++i)
        if (xs[i] < -M_PI/2) xs[i] += 2*M_PI;
    }
  }
  calc_ext();
  if (max_lon > 2*M_PI - min_lon)
    for (int i = 0; i < ne; ++i) xs[i] -= 2*M_PI;
  if (-min_lon > 2*M_PI + max_lon)
    for (int i = 0; i < ne; ++i) xs[i] += 2*M_PI;  
}

// [[[xs]*nvmax [ys]*nvmax]*]
extern "C" void get_oriented_cell_vertex_lonlat_list (Real* xys, int wrap) {
  assert(s_mm);
  const auto m = *s_mm;
  ompparfor for (int ci = 0; ci < m.t.ncell; ++ci) {
    Real* const xs = xys + 2*m.t.max_vtx*ci;
    Real* const ys = xs + m.t.max_vtx;
    int elem[max_vtx];
    get_oriented_cell(m.t, m.v, ci, elem);
    const int ne = m.t.c2ne[ci];
    for (int i = 0; i < ne; ++i) {
      const int ei = elem[i];
      assert(ei >= 0 && ei < m.t.nvtx);
      xyz2lonlat(m.v.x[ei], m.v.y[ei], m.v.z[ei], xs[i], ys[i]);
    }
    wrap_lonlat(xs, ys, ne);
  }
}

// [[[xs]*nvmax [ys]*nvmax]*]
extern "C" void
get_oriented_cell_vertex_xy_list (const Real X, const Real Y, Real* xys) {
  assert(s_mm);
  const auto m = *s_mm;
  ompparfor for (int ci = 0; ci < m.t.ncell; ++ci) {
    Real* const xs = xys + 2*m.t.max_vtx*ci;
    Real* const ys = xs + m.t.max_vtx;
    int elem[max_vtx];
    get_oriented_cell(m.t, m.v, ci, elem);
    const int ne = m.t.c2ne[ci];
    Real xmin = 1.1*X, xmax = 0, ymin = 1.1*Y, ymax = 0;
    for (int i = 0; i < ne; ++i) {
      const int ei = elem[i];
      assert(ei >= 0 && ei < m.t.nvtx);
      xs[i] = m.v.x[ei];
      ys[i] = m.v.y[ei];
      xmin = std::min(xmin, xs[i]); xmax = std::max(xmax, xs[i]);
      ymin = std::min(ymin, ys[i]); ymax = std::max(ymax, ys[i]);
    }
    if (xmax - xmin > 0.5*X)
      for (int i = 0; i < ne; ++i)
        if (xs[i] < 0.5*X) xs[i] += X;
    if (ymax - ymin > 0.5*Y)
      for (int i = 0; i < ne; ++i)
        if (ys[i] < 0.5*Y) ys[i] += Y;
  }
}

extern "C" void make_yz_image (const int nz /* fast index */, const Real* f3d,
                               int* ony, Real* f2d) {
  assert(s_mm);
  const auto m = *s_mm;

  std::vector<Real> suy; {
    std::set<Real> uy;
    for (int i = 0; i < m.t.ncell; ++i) uy.insert(m.c.y[i]);
    for (const auto e : uy) suy.push_back(e);
  }
  *ony = suy.size();
  const int ny = *ony;

  for (int i = 0; i < ny*nz; ++i) f2d[i] = 0;
  std::vector<int> cnt(ny, 0);
  for (int iy = 0; iy < m.t.ncell; ++iy) {
    const int siy = (std::lower_bound(suy.begin(), suy.end(), m.c.y[iy])
                     - suy.begin());
    assert(siy >= 0 && siy < ny);
    ++cnt[siy];
    for (int iz = 0; iz < nz; ++iz)
      f2d[nz*siy + iz] += f3d[nz*iy + iz];
  }
  for (int iy = 0; iy < ny; ++iy)
    for (int iz = 0; iz < nz; ++iz)
      f2d[nz*iy + iz] /= cnt[iy];
}

extern "C" void advect_ir (int output, const char* optstr, int nstep,
                           Real* ic, Real* mid, Real* fin) {
  assert(s_mm);
  ir::TestProblemOptions o;
  o.flow_ic_str = optstr;
  o.nstep = nstep;
  o.movie = output == 2;
  const bool all = output == 1;
  const auto he = Comm::create();
  ir::Params p;
  ir::set(*s_mm, *he, o, p);
  const auto op = ir::RemapFormOp::create(s_mm, s_m, he, p);
  ir::run_test_problem(op, o, nullptr, (all || o.movie) ? ic : nullptr,
                       all ? mid : nullptr, all ? fin : nullptr);
}

extern "C" void advect_mesh (const char* optstr, const int nstep, const int stepno,
                             Real* lls_orig, Real* lls_advected) {
  assert(s_mm);
  ir::TestProblemOptions o;
  o.flow_ic_str = optstr;
  o.nstep = nstep;
  o.movie = false;
  std::vector<Real> xyzs_orig, xyzs_advected;
  ir::advect_mesh(s_mm, s_m, o, stepno, xyzs_orig, xyzs_advected);
  const auto& m = *s_mm;
  const int
    fac = optstr[2] == '3' ? 2 : 1,
    nper1d = fac*m.t.max_vtx, nper2d = 2*nper1d, nper3d = 3*nper1d;
  for (int ci = 0; ci < m.t.ncell; ++ci) {
    for (int mesh = 0; mesh < 2; ++mesh) {
      const Real
        * base_orig = mesh == 0 ? xyzs_orig.data() : xyzs_advected.data(),
        * xs = base_orig + ci*nper3d,
        * ys = xs + nper1d,
        * zs = ys + nper1d;
      Real
        * base_advected = mesh == 0 ? lls_orig : lls_advected,
        * lon = base_advected + ci*nper2d,
        * lat = lon + nper1d;
      const int ne = fac*m.t.c2ne[ci];
      for (int i = 0; i < ne; ++i)
        xyz2lonlat(xs[i], ys[i], zs[i], lon[i], lat[i]);
      wrap_lonlat(lon, lat, ne);
    }
  }
}

extern "C" void unittest_ir (int output, Real* data) {
  assert(s_mm);
  const auto he = Comm::create();
  ir::unittest(s_mm, s_m, he, output ? data : nullptr);
}

extern "C" void make_graph_file (const char* fname) {
  assert(s_mm);
  const auto& m = *s_mm;
  const auto fp = FilePtr(fopen(fname, "w"));
  fprintf(fp.get(), "%d %d\n", m.t.ncell, m.t.nedge);
  for (int ci = 0; ci < m.t.ncell; ++ci) {
    for (int j = 0; j < m.t.c2ne[ci]; ++j) {
      const auto ej = m.t.get_c2e(ci,j);
      const auto cj = get_other_cell(m.t, ci, ej);
      fprintf(fp.get(), " %d", cj+1);
    }
    fprintf(fp.get(), "\n");
  }
}

extern "C" void vert_lat_slice_get_cols (
  const PyInt ncgid,
  const PyInt* const cgids, // cell GIDs in the order that field data are ordered
  Real* const lons_list)    // pairs of (min-lon, max-lon) extents of a cell
{
  assert(s_mm);
  const auto& m = *s_mm;
  const auto& v = m.v;
  for (int k = 0; k < ncgid; ++k) {
    const Int ci = cgids[k];
    const int ne = m.t.c2ne[ci];
    Real lon[max_vtx], lat[max_vtx];
    for (int i = 0; i < ne; ++i) {
      const int vi = m.t.get_c2v(ci,i);
      xyz2lonlat(v.x[vi], v.y[vi], v.z[vi], lon[i], lat[i]);
    }
    wrap_lonlat(lon, lat, ne);
    Real lon_min = 10, lon_max = -10;
    for (int i = 0; i < ne; ++i) {
      lon_min = std::min(lon_min, lon[i]);
      lon_max = std::max(lon_max, lon[i]);
    }
    lons_list[2*k  ] = lon_min;
    lons_list[2*k+1] = lon_max;
  }
}

// For a plane given by axis, develop an angular coordinate theta and provide
// theta extents for each cell in the cgids list.
extern "C" void vert_axis_slice_get_cols (
  const PyInt ncgid, const Vec axis_,
  const PyInt* const cgids, // cell GIDs in the order that field data are ordered
  Real* const thetas_list)  // pairs of (min-theta, max-theta) extents of a cell
{
  assert(s_mm);
  const auto& m = *s_mm;
  const auto& v = m.v;
  Vec axis; geo::copy(axis, axis_); geo::normalize(axis);
  Vec xaxis, yaxis, ref = {0};
  if (std::abs(axis[2]) < 1/std::sqrt(2)) ref[2] = 1;
  else ref[0] = 1;
  geo::cross(ref, axis, yaxis); geo::normalize(yaxis);
  geo::cross(axis, yaxis, xaxis); geo::normalize(xaxis);
  for (int k = 0; k < ncgid; ++k) {
    const Int ci = cgids[k];
    assert(ci >= 0);
    assert(ci < m.t.ncell);
    const int ne = m.t.c2ne[ci];
    assert(ne <= max_vtx);
    Real thetas[max_vtx], theta_min = 10, theta_max = -10;
    for (int i = 0; i < ne; ++i) {
      const int vi = m.t.get_c2v(ci,i);
      assert(vi >= 0);
      assert(vi < m.t.nvtx);
      Vec pos = {v.x[vi], v.y[vi], v.z[vi]};
      thetas[i] = std::atan2(geo::dot(pos, yaxis), geo::dot(pos, xaxis));
      theta_min = std::min(theta_min, thetas[i]);
      theta_max = std::max(theta_max, thetas[i]);
    }
    if (theta_max - theta_min > M_PI)
      for (int i = 0; i < ne; ++i)
        if (thetas[i] > 0)
          thetas[i] -= 2*M_PI;
    theta_min = 10; theta_max = -10;
    for (int i = 0; i < ne; ++i) {
      theta_min = std::min(theta_min, thetas[i]);
      theta_max = std::max(theta_max, thetas[i]);
    }
    thetas_list[2*k  ] = theta_min;
    thetas_list[2*k+1] = theta_max;
  }
}

namespace {
struct OdeFn : public Trajectory::OdeFn {
  mutable std::vector<Real> ts, lons, lats, zetas;
  OdeFn (const VelocityField& vf_) : Trajectory::OdeFn(vf_, 0) {}
  void record (Real t, const Real* y) const override {
    ts.push_back(t); lons.push_back(y[0]); lats.push_back(y[1]);
    if (zcomp) zetas.push_back(y[2]);
  }
};
} // namespace

extern "C" void integrate_gallery_velocity_field (
  const char vf_short, const Int n_space_pts, const Int n_time_pts,
  Real* ts, // Record time.
  // 2D arrays. Fast index is time points. First time point contains ICs.
  Real* lons, Real* lats, Real* zetas)
{
  if (vf_short == 'a') {
    printf("vf_short 'a' is not supported\n");
    return;
  }
  const auto vf = vf_from_short(vf_short);
  if ( ! vf) {
    printf("vf_short %c is invalid\n", vf_short);
    return;
  }
  const bool zcomp = vf->is_3d();
  ompparfor for (int si = 0; si < n_space_pts; ++si) {
    timeint::Options opts;
    opts.set_abs_tol(std::numeric_limits<Real>::epsilon());
    opts.set_rel_tol(1e-10);
    Real pos[3], wrk[27];
    pos[0] = lons[n_time_pts*si];
    pos[1] = lats[n_time_pts*si];
    if (zcomp) pos[2] = zetas[n_time_pts*si];
    timeint::Info info;
    OdeFn f(*vf);
    timeint::ark45(opts, f, pos, zcomp ? 3 : 2, 0, vf->cycle_time(), wrk, &info);
    const int n = static_cast<int>(f.ts.size());
    const int os = n_time_pts*si;
    const auto copy = [&] (const int d, const int s) {
      assert(os+d >= 0); assert(os+d < n_space_pts*n_time_pts);
      assert(s >= 0); assert(s < n);
      ts[os+d] = f.ts[s];
      lons[os+d] = f.lons[s];
      lats[os+d] = f.lats[s];
      if (f.zcomp) zetas[os+d] = f.zetas[s];
    };
    Real progress = 0;
    int k = 0;
    for (int i = 0; i < n; ++i) {
      if (Real(i)/(n-1) < progress) continue;
      copy(k, i);
      if (++k == n_time_pts) break;
      progress = Real(k)/(n_time_pts-1);
    }
    if (k < n_time_pts) {
      for ( ; k < n_time_pts; ++k) copy(k, n-1);
    }
    // Be sure last time point is included.
    copy(n_time_pts-1, n-1);
  }
}
