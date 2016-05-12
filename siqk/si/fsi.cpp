#include "Array_raw.hpp"
#include "sik.hpp"

extern "C" void clipagainstpolysphere_ (
  double const* const clip_poly, int const* const clip_poly_n_vertices,
  double const* const clip_edge_normals, double const* const vi, int const* const ni,
  double* const vo, int* const no, double* const wrk, int const* const n_vertices,
  int* const info)
{
  Array2D<double> avo(3, *n_vertices, vo);
  const bool success = siqk::sh::clip_against_poly<siqk::SphereGeometry>(
    Array2D<const double>(3, *clip_poly_n_vertices, clip_poly),
    Array2D<const double>(3, *clip_poly_n_vertices, clip_edge_normals),
    Array2D<const double>(3, *ni, vi), *ni,
    avo, *no, wrk, *n_vertices);
  *info = success ? 0 : 1;
}
