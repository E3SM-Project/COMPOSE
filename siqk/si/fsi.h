// Fortran interface to simple polygon clipping routine.

extern "C" void clipagainstpolysphere_(
  // 3 x clip_poly_n_vertices clip spherical polygon vertex list.
  double const* const clip_poly, int const* const clip_poly_n_vertices,
  // 3 x clip_poly_n_vertices clip polygon's inward-facing edge normals.
  double const* const clip_edge_normals,
  // 3 x ni polygon to clip.
  double const* const to_clip_poly, int const* const ni,
  // On output, a 3 x no clipped polygon.
  double* const vo, int* const no,
  // Workspace. Both vo and wrk must have n_vertices of space available.
  double* const wrk, int const* const n_vertices,
  // info = 0 on success. info = 1 if n_vertices is not large enough.
  int* const info);
