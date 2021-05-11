#include <typeinfo>

#include "slmm_c_compat.hpp"

#include "slmm_spf.hpp"
#include "slmm_gll.hpp"
#include "slmm_mesh.hpp"
#include "slmm_util.hpp"

void init () {
  static bool inited = false;
  if ( ! inited) Kokkos::initialize();
  inited = true;
}

void finalize () {
  static bool fin = false;
  if ( ! fin) Kokkos::finalize();
  fin = true;
}

int solve_1eq_bc_qp_c (
  const int n, const double* w, const double* a, const double b,
  const double* xlo, const double* xhi, const double* y, double* x,
  const int max_its)
{
  return slmm::spf::solve_1eq_bc_qp(
    n, w, a, b, xlo, xhi, false, y, x, max_its);
}

void make_mesh (const int ne, const int np, const int nonuni,
                double* const p, int* const e) {
  slmm::AVec3s geo_p;
  slmm::AIdxs geo_c2n;
  slmm::mesh::make_cubedsphere_mesh(geo_p, geo_c2n, ne);
  if (nonuni) slmm::mesh::make_nonuniform(geo_p);
  if (np > 0) {
    slmm::AVec3s cgll_p;
    slmm::AIdxs cgll_c2n, cgll_io_c2n;
    slmm::mesh::make_cgll_from_geo(geo_p, geo_c2n, np, slmm::GLL(), cgll_p, cgll_c2n);
    slmm::mesh::make_io_cgll_from_internal_cgll(cgll_p, cgll_c2n, cgll_io_c2n);
    std::copy(cgll_p.data(), cgll_p.data() + cgll_p.size(), p);
    std::copy(cgll_io_c2n.data(), cgll_io_c2n.data() + cgll_io_c2n.size(), e);
  } else {
    std::copy(geo_p.data(), geo_p.data() + geo_p.size(), p);
    std::copy(geo_c2n.data(), geo_c2n.data() + geo_c2n.size(), e);
  }
}

static void wrap_lonlat (double* xs, double* ys, const int ne) {
  double min_lon, max_lon, min_lat, max_lat;
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

void make_latlon_quads (const int ncell, const double* ps, const int* es, double* lls) {
  static constexpr int nv = 4;
  for (int ie = 0; ie < ncell; ++ie) {
    const auto e = es + nv*ie;
    const auto lons = lls + 2*nv*ie;
    const auto lats = lons + 4;
    for (int iv = 0; iv < nv; ++iv) {
      const auto p = ps + 3*e[iv];
      slmm::xyz2ll(p[0], p[1], p[2], lats[iv], lons[iv]);
    }
    wrap_lonlat(lons, lats, nv);
  }
}
