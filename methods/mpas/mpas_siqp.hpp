// Derived from SIQK package of https://github.com/E3SM-Project/COMPOSE
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef MPAS_INCLUDE_SIQ_HPP
#define MPAS_INCLUDE_SIQ_HPP

#include <cassert>
#include <cmath>
#include <memory>
#include <type_traits>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <cstring>

#ifdef SIQP_TIME
# include <unistd.h>
# include <sys/time.h>
# include <sys/resource.h>
#endif

#include "mpas_array.hpp"

namespace siqp {
typedef double Real;

#ifdef SIQP_TIME
timeval tic();
double calc_et(const timeval& t1, const timeval& t2);
double toc(const timeval& t1);
double get_memusage();
#else
inline int tic () { return 0; }
inline double toc (const int&) { return 0; }
#endif
void print_times(const std::string& name, const double* const parts,
                 const int nparts);

template <typename V, typename CV>
void copy (V dst, CV src, const int n) {
  for (int i = 0; i < n; ++i) dst[i] = src[i];
}

// A decorator function so that a for loop's counter can be auto typed.
template <typename V>
typename V::size_type zero (const V& v) { return 0; }

template <typename CV, typename V>
int solve_quadratic_equation (const CV a, const CV b, const CV c, V xs) {
  const double disc = b*b - 4*a*c;
  if (disc < 0) return 0;
  if (disc == 0) {
    xs[0] = -0.5*b/a;
    return 1;
  }
  if (b >= 0) {
    const double t = b + std::sqrt(disc);
    xs[0] = -0.5*t/a;
    xs[1] = -2*c/t;
  } else {
    const double t = std::sqrt(disc) - b;
    xs[0] = 2*c/t;
    xs[1] = 0.5*t/a;
  }
  return 2;
}

inline void error (const std::string& msg)
{ throw std::runtime_error(msg.c_str()); }

// Sutherland-Hodgmann polygon clipping algorithm. Follow Foley, van Dam,
// Feiner, Hughes Fig 3.49.
namespace sh {
/* A mesh is described by the following arrays:
       p: 3 x #nodes, the array of vertices.
       e: max(#verts) x #elems, the array of element base-0 indices.
       nml: 3 x #edges, the array of edge normals.
       en: max(#verts) x #elems, the array of edge-normal base-0 indices.
     e. e indexes p. e(i,j) == -1 in column j indicates that j:end are not used.
     nml. As a mesh is refined, cancellation error makes an edge normal based
   off of an element's vertices increasingly inaccurate. Roughly, if an edge
   subtends angle phi of the sphere, -log10(phi/(2 pi)) digits are lost in the
   edge normal. Therefore, we compute edge normals offline, since in certain
   meshes, they can be computed by an accurate means. E.g., in a cubed-sphere
   mesh, the whole line of a square face can be used to compute the edge
   normal. Furthermore, there are far fewer unique edge normals than edges.
 */
struct Mesh {
  typedef std::shared_ptr<Mesh> Ptr;
  typedef std::shared_ptr<const Mesh> ConstPtr;
  Array2D<const Real> p, nml;
  Array2D<const int> e, en;
};

// Efficient user routine that uses the mesh data structure.
bool clip_against_poly(
  // Clip mesh. m.e(:,cp_e) is the element, and m.en(:,cp_e) is the
  // corresponding list of normal indices.
  const Mesh& m, const int cp_e,
  // A list of vertices describing the polygon to clip. The vertices must be in
  // a convention-determined order, such as CCW. vi(:,1:ni-1) are valid entries.
  const Array2D<const Real>& vi, const int ni,
  // On output, vo(:,0:no-1) are vertices of the clipped polygon. no is 0 if
  // there is no intersection.
  Array2D<Real>& vo, int& no,
  // Workspace. nvertwrk applies to both wrk and vo.n(). If nvertwrk is not
  // large enough, false is returned. wrk should have at least 3*nvertwrk slots.
  Real* const wrk, const int nvertwrk);

// ICE: Intersection with curved edges.
//
// Some terminology:
//   s, m, p: start, middle, end points of a curved edge. m is not really the
//     middle or midpoint; indeed, it is unlikely to be on the curve. Rather,
//     it's a point that defines the curve.
//   smp-plane: plane defined by s, m, p.
//   straight: straight on a plane, or a great arc on the sphere.
//   curved: quadratic on a plane, projected quadratic on the sphere.
//   cedge: a curved edge.
//   sedge: a straight edge, including a great arc.
//   cpoly, spoly: similar terminology; but note that a cpoly can contain a mix
//     of cedges and sedges.
//   ed: short for edge
//   vt: node type. s, p are endpoint vertices (0); m is a midpoint node (1).
//   ced is [s m p].
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
//     x(1/2) = M => m = 4 M - s - p,
// where M is a point that is intended to be on the curve and serves as a useful
// midpoint reference. This construction has the essential property that the
// curve is invariant to the swapping of s and p. However, the clip routines are
// independent of this definition; (s,m,p) are used as in equation (1), and that
// is all that is needed.
//   When segments are extracted from x(a) in a clip, we use c in the segment
// [0,c], d in the segment [d,1], and both in [c,d]. Similarly, x(c) = r is the
// new start point, and x(d) = q is the new end point. A segment requires a
// midpoint so that the resulting parameterized curve sits on the original; this
// is n. Hence a segment of (s,m,p) defined by [c,d] subset [0,1] is (r,n,q). n
// is given by
//     n = 2 (c d - c - d + 1) s + (c + d - 2 c d) m + 2 c d p,
// which satisfies
//     x(a) = (1-a)^2 s + a(1-a) m + a^2 p = (1-b)^2 r + b(1-b) n + b^2 q
// for all b in [0,1], where b = (a-c)/(d-c), r = x(c), q = x(d).
bool clip_cpoly_against_convex_poly(
  // Clip mesh. m.e(:,cp_e) is the element, and m.en(:,cp_e) is the
  // corresponding list of normal indices.
  const Mesh& m, const int cp_e,
  // Input vertex and vertex type lists.
  const Array2D<const Real>& vi, const Array1D<const int>& vti, const int ni,
  // Outputs.
  Array2D<Real>& vo, Array1D<int>& vto, int& no,
  // Workspace. nwrk applies to both *wrk and v*o.  If workspace is not large
  // enough, false is returned. For real quantities, size is 3*nwrk; for
  // integer, nwrk.
  Real* const rwrk, int* const iwrk, const int nwrk);

void sedge_eval(const Real s[3], const Real p[3], const Real& a, Real u[3]);
void cedge_eval(const Real s[3], const Real m[3], const Real p[3], const Real& a,
                Real u[3]);

// Determine m in (s,m,p) such that n is on the cuve. n and m may be aliases.
void middle_matches(const Real s[3], const Real p[3], const Real n[3], Real m[3]);

// A vctri is a triangle with 2 straight edges and one curve. The curved is the
// middle edge, so vts = [0 0 1 0]. s is v2 and p is v3.
Real calc_vctri_jacobian(const Real v1[3], const Real s[3], const Real m[3], const Real p[3],
                         const Real alpha[3], Real u[3], const bool nosign = true);

bool test_ice_math();
} // namespace sh

template <typename T> constexpr T square (const T& x) { return x*x; }

// Octree for search.
class Octree {
public:
  typedef Real BoundingBox[6];
  typedef std::shared_ptr<Octree> Ptr;
  typedef std::shared_ptr<const Octree> ConstPtr;

  struct Options {
    // Do not go beyond max_depth depth, including the root and leaf. With this
    // constraint, try to go deep enough so that a leaf has no more than
    // max_nelem elements.
    int max_depth, max_nelem;
    Options () : max_depth(10), max_nelem(8) {}
  };

  struct Functor {
    virtual ~Functor () {}
    virtual void operator() (const int element_index) = 0;
  };

  // Bounding box for a cluster of points ps (possibly vertices).
  template <typename Ar2> static void
  calc_bb(const Ar2& ps, const int np, BoundingBox bb);
  template <typename Ar2> static void
  calc_bb(const Ar2& ps, BoundingBox bb);
  template <typename Ar2> static void
  calc_bb(const Ar2& p, const Array2D<const int>& e, Array2D<Real>& ebbs);
  template <typename Ar2> static void
  calc_bb(const Ar2& p, const int* e, const int ne, Real* ebb);

  // p is a 3xNp array of points. e is a KxNe array of elements. An entry <0 is
  // ignored. All <0 entries must be at the end of an element's list.
  Octree(const Array2D<const Real>& p, const Array2D<const int>& e,
         const Options& o);
  Octree(const Array2D<const Real>& p, const Array2D<const int>& e);

  // Apply f to every element in leaf nodes with which bb overlaps.
  void apply(const BoundingBox& bb, Functor& f) const;

private:
  /* Each node in the oct-tree contains 8 integers, stored in 'nodes'.

     >0 is an index into 'nodes', pointing to a child node.

     A <=0 entry in 'nodes' indicates a leaf node. If 0, there are no elements
     in the leaf. If <0, the negative of the entry minus 1 is the index of an
     offset array indexing 'elems'.

     Each segment of 'elems' contains a list of element indices covered by a
     leaf node. Element indices refer to the list of elements the caller
     provides during oct-tree construction.
  */

  // nodes(:,i) is a list. The list includes children of node i (>0) and leaf
  // node data (<=0).
  Array2D<int> nodes_;
  // A leaf node corresponding to -k covers elements
  //     elems[offset[k] : offset[k]-1].
  Array1D<int> offset_, elems_;
  // Root node's bounding box.
  BoundingBox bb_;

  void init(const Array2D<const Real>& p, const Array2D<const int>& e,
            const Options& o);
  void apply_r(const int ni, const BoundingBox& nbb, const BoundingBox& bb,
               Functor& f) const;
};

// If a bounding box was constructed from vertices of a spherical polygon,
// expand it to account for the possible protrusion of the sphere.
template <typename BB>
static void pad_bb (BB bb) {
  Real hl = 0.5*std::sqrt(square(bb[3] - bb[0]) + square(bb[4] - bb[1]) +
                          square(bb[5] - bb[2]));
  // Limit the half-length to the circle's radius.
  hl = std::min(1.0, hl);
  // Max distance from a chord of length 2 hl to the unit circle:
  //     hl = sin theta
  //    pad = 1 - cos theta = 1 - sqrt(1 - sin^2 theta) = 1 - sqrt(1 - hl^2).
  const Real pad = 1 - std::sqrt(1 - square(hl));
  for (int i = 0; i < 3; ++i) bb[  i] -= pad;
  for (int i = 0; i < 3; ++i) bb[3+i] += pad;
}

template <typename Ar2> void Octree
::calc_bb (const Ar2& p, const int* e, const int ne, Real* ebb) {
  for (int j = 0; j < 3; ++j)
    ebb[j] = ebb[j+3] = p(j, e[0]);
  for (int i = 1; i < ne; ++i) {
    if (e[i] == -1) break;
    for (int j = 0; j < 3; ++j) {
      ebb[j] = std::min(ebb[j], p(j, e[i]));
      ebb[j+3] = std::max(ebb[j+3], p(j, e[i]));
    }
  }
  pad_bb(ebb);
}

template <typename Ar2> void Octree
::calc_bb (const Ar2& ps, const int np, BoundingBox bb) {
  if (np == 0) return;
  for (int j = 0; j < 3; ++j)
    bb[j] = bb[j+3] = ps(j,0);
  for (int i = 1; i < np; ++i)
    for (int j = 0; j < 3; ++j) {
      bb[j] = std::min(bb[j], ps(j,i));
      bb[j+3] = std::max(bb[j+3], ps(j,i));
    }
  pad_bb(bb);
}

template <typename Ar2> void Octree::calc_bb (const Ar2& ps, BoundingBox bb) {
  calc_bb(ps, ps.n(), bb);
}

template <typename Ar2> void Octree
::calc_bb (const Ar2& p, const Array2D<const int>& e, Array2D<Real>& ebbs) {
  assert(ebbs.n() == e.n());
  for (int k = 0; k < e.n(); ++k)
    calc_bb(p, e(k), e.m(), ebbs(k));
}

// Fill en and nml arrays. Normals are inward.
void fill_normals(sh::Mesh& m);

namespace test {
template <typename Functor, typename Scalar>
void parallel_reduce (const int n, Functor f, Scalar& r) {
  for (int i = 0; i < n; ++i)
    f(i, r);
}

// Calculates the area given a polygon with n vertices v.
struct AreaFunctor {
  virtual ~AreaFunctor () {}
  virtual Real operator() (const Real* const v, const int n) const = 0;
};

Real test_area_ot (const Array2D<const Real>& cp, const Array2D<const int>& ce,
                   const Array2D<const Real>& p, const Array2D<const int>& e,
                   const bool verbose, const bool use_ice,
                   // Optional. If one is not provide, use a bad formula.
                   AreaFunctor* af = nullptr);
} // namespace test
} // namespace siqp

#endif
