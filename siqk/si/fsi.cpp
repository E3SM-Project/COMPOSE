#include "Array_raw.hpp"
#include "siqp.hpp"

extern "C" void clipagainstpolysphere_ (
  double const* const clip_poly, int const* const clip_poly_n_vertices,
  double const* const clip_edge_normals, double const* const vi, int const* const ni,
  double* const vo, int* const no, double* const wrk, int const* const n_vertices,
  int* const info)
{
  Array2D<double> avo(3, *n_vertices, vo);
  const bool success = siqp::sh::clip_against_poly<siqp::SphereGeometry>(
    Array2D<const double>(3, *clip_poly_n_vertices, clip_poly),
    Array2D<const double>(3, *clip_poly_n_vertices, clip_edge_normals),
    Array2D<const double>(3, *ni, vi), *ni,
    avo, *no, wrk, *n_vertices);
  *info = success ? 0 : 1;
}

extern "C" void iceclipagainstpolyplane_(
  double const* const clip_poly, int const* const clip_poly_n_vertices,
  double const* const clip_edge_normals,
  double const* const vi, int const* const vti, int const* const ni,
  double* const vo, int* const vto, int* const no,
  double* const rwrk, int* const iwrk, int const* const nwrk,
  int* const info)
{
  Array2D<double> avo(3, *nwrk, vo);
  Array1D<int> avto(*nwrk, vto);
  const bool success =
    siqp::ice<siqp::PlaneGeometry>::clip_cpoly_against_convex_poly(
      Array2D<const double>(3, *clip_poly_n_vertices, clip_poly),
      Array2D<const double>(3, *clip_poly_n_vertices, clip_edge_normals),
      Array2D<const double>(3, *ni, vi), Array1D<const int>(*ni, vti), *ni,
      avo, avto, *no, rwrk, iwrk, *nwrk);
  *info = success ? 0 : 1;
}
