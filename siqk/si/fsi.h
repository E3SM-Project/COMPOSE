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

// ICE: Intersection with curved edges.
//
// Some terminology:
//   s, m, p: start, middle, end points of a curved edge. m is not really the
//     middle or midpoint; indeed, it is unlikely to be on the curve. Rather,
//     it's a point that defines the curve. See below for more.
//   straight: straight on a plane, or a great arc on the sphere.
//   curved: quadratic on a plane, projected quadratic on the sphere.
//   cedge: a curved edge.
//   sedge: a straight edge, including a great arc.
//   cpoly, spoly: similar terminology; but note that a cpoly can contain a mix
//     of cedges and sedges.
//   (vs, vts, n): Vertex list. vs is an array of vertices. There are n
//     vertices. vts is a list of vertex types. s, p are endpoint vertices (0);
//     m is a midpoint node (1). If an edge is straight, then vts(k:k+1) = [0
//     0]; if an edge is curved, then vts(k:k+2) = [0 1 0]. Keep in mind that
//     there is an edge that wraps around the end of the list. The wrap can
//     occur like 0|1 0 or like 0 1|0. As an example, [1 0 0 1 0] is a vertex
//     type list for a triangle containing two cedges and one sedge.
//
// Some math.
//   a in [0,1] is the parameter in the curve
//     x(a) = (1-a)^2 s + a (1-a) m + a^2 p.                                 (1)
//   s and p sit on the curve, but m in general does not. We can define m by
//     x(1/2) = M => m = 4 n - s - p,
// where M is a point that is intended to be on the curve and serves as a useful
// midpoint reference. This construction has the essential property that the
// curve is invariant to the swapping of s and p. However, the clip routines are
// independent of this definition; (s,m,p) are used as in equation (1), and that
// is all that is needed.
//   When segments are extracted from x(a) in a clip, we use c in the segment
// [c,1], d in the segment [0,d], and both in [c,d]. Similarly, x(c) = r is the
// new start point, and x(d) = q is the new end point. A segment requires a
// midpoint so that the resulting parameterized curve sits on the original; this
// is n. Hence a segment of (s,m,p) defined by [c,d] subset [0,1] is (r,n,q). n
// is given by
//     n = 2 (c d - c - d + 1) s + (c + d - 2 c d) m + 2 c d p,
// which satisfies
//     x(a) = (1-a)^2 s + a(1-a) m + a^2 p = (1-b)^2 r + b(1-b) n + b^2 q
// for all b in [0,1], where b = (a-c)/(d-c), r = x(c), q = x(d).
//
// Even though this is a 2D routine, the vertices and vectors must still be
// 3D. The third value is ignored.
extern "C" void iceclipagainstpolyplane_(
  // 3 x clip_poly_n_vertices clip spherical polygon vertex list.
  double const* const clip_poly, int const* const clip_poly_n_vertices,
  // 3 x clip_poly_n_vertices clip polygon's inward-facing edge normals.
  double const* const clip_edge_normals,
  // 3 x ni curved polygon to clip.
  double const* const vi, int const* const vti, int const* const ni,
  // On output, a 3 x no clipped polygon.
  double* const vo, int* const vto, int* const no,
  // Workspace. rwrk is 3*nwrk; iwrk is 1*nwrk.
  double* const rwrk, int* const iwrk, int const* const nwrk,
  // info = 0 on success. info = 1 if workspace is not large enough.
  int* const info);
