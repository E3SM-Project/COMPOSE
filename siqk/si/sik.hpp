#ifndef INCLUDE_SIK_HPP
#define INCLUDE_SIK_HPP

#include <cmath>
#include <cassert>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <set>

#ifdef SIQK_TIME
# include <unistd.h>
# include <sys/time.h>
# include <sys/resource.h>
#endif

namespace siqk {
#ifdef SIQK_TIME
static timeval tic () {
  timeval t;
  gettimeofday(&t, 0);
  return t;
}
static double calc_et (const timeval& t1, const timeval& t2) {
  static const double us = 1.0e6;
  return (t2.tv_sec * us + t2.tv_usec - t1.tv_sec * us - t1.tv_usec) / us;
}
static double toc (const timeval& t1) {
#ifdef SIQK_USE_KOKKOS
  Kokkos::fence();
#endif
  timeval t;
  gettimeofday(&t, 0);
  return calc_et(t1, t);
}
static double get_memusage () {
  static const double scale = 1.0 / (1 << 10); // Memory in MB.
  rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  return ru.ru_maxrss*scale;
}
#else
static inline int tic () { return 0; }
static inline double toc (const int&) { return 0; }
#endif
static void print_times (const std::string& name, const double* const parts,
                         const int nparts) {
#ifdef SIQK_TIME
  double total = 0; for (int i = 0; i < nparts; ++i) total += parts[i];
  printf("%20s %1.3e s %7.1f MB", name.c_str(), total, get_memusage());
  for (int i = 0; i < nparts; ++i) printf(" %1.3e s", parts[i]);
  printf("\n");
#endif
}

template <typename V, typename CV>
static void copy (V dst, CV src, const int n) {
  for (int i = 0; i < n; ++i) dst[i] = src[i];
}

// A decorator function so that a for loop's counter can be auto typed.
template <typename V> KOKKOS_INLINE_FUNCTION
typename V::size_type zero(const V& v) { return 0; }

// Planar geometry calculations.
struct PlaneGeometry {
  template <typename V> KOKKOS_INLINE_FUNCTION
  static void scale (const double a, V v) {
    v[0] *= a; v[1] *= a;
  }
  template <typename CV> KOKKOS_INLINE_FUNCTION
  static double dot_c_amb (const CV c, const CV a, const CV b) {
    return c[0]*(a[0] - b[0]) + c[1]*(a[1] - b[1]);
  }
  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void combine (const CV u, const CV v, const double a, V x) {
    const double oma = 1 - a;
    x[0] = oma*u[0] + a*v[0];
    x[1] = oma*u[1] + a*v[1];
  }

  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void edge_normal (const CV e1, const CV e2, V en) {
    en[0] = e1[1] - e2[1];
    en[1] = e2[0] - e1[0];
  }

  template <typename CV> KOKKOS_INLINE_FUNCTION
  static bool inside (const CV v, const CV e1, const CV e2, const CV en) {
    return dot_c_amb(en, v, e1) > 0 && dot_c_amb(en, v, e2) > 0;
  }

  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void intersect (const CV v1, const CV v2, const CV e1, const CV en,
                         V intersection) {
    double a; {
      const double
        num = dot_c_amb(en, e1, v1),
        den = dot_c_amb(en, v2, v1);
      a = num == 0 || den == 0 ? 0 : num/den;
      a = a < 0 ? 0 : a > 1 ? 1 : a;
    }
    combine(v1, v2, a, intersection);
  }

  template <typename CV> KOKKOS_INLINE_FUNCTION
  static bool output (const CV v, int& no, Array2D<double>& vo) {
#ifdef SIKQ_DEBUG
    if (no >= vo.n()) {
      std::stringstream ss;
      ss << "output: No room in vo; vo.n() is " << vo.n() << " but no is "
         << no << "\n";
      error(ss.str().c_str());
    }
#endif
    if (no >= vo.n()) return false;
    vo(0,no) = v[0];
    vo(1,no) = v[1];
    ++no;
    return true;
  }

  //todo Handle non-convex case.
  KOKKOS_INLINE_FUNCTION
  static double calc_area (const Array2D<const double>& v) {
    double area = 0;
    for (int i = 1; i < v.n() - 1; ++i) {
      double v1[2], v2[2];
      v1[0] = v(0,i) - v(0,0);
      v1[1] = v(1,i) - v(1,0);
      v2[0] = v(0,i+1) - v(0,0);
      v2[1] = v(1,i+1) - v(1,0);
      const double a = v1[0]*v2[1] - v1[1]*v2[0];
      area += a;
    }
    return 0.5*area;
  }
};

// Geometry on the sphere. All inputs and outputs are relative to the
// unit-radius sphere.
struct SphereGeometry {
  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void cross (const CV a, const CV b, V c) {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
  }
  template <typename CV> KOKKOS_INLINE_FUNCTION
  static double dot (const CV a, const CV b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  }
  template <typename CV> KOKKOS_INLINE_FUNCTION
  static double norm2 (const CV v) {
    return dot(v, v);
  }
  template <typename V> KOKKOS_INLINE_FUNCTION
  static void scale (const double a, V v) {
    v[0] *= a; v[1] *= a; v[2] *= a;
  }
  template <typename V> KOKKOS_INLINE_FUNCTION
  static void normalize (V v) {
    scale(1.0/std::sqrt(norm2(v)), v);
  }
  template <typename CV> KOKKOS_INLINE_FUNCTION
  static double dot_c_amb (const CV c, const CV a, const CV b) {
    return c[0]*(a[0] - b[0]) + c[1]*(a[1] - b[1]) + c[2]*(a[2] - b[2]);
  }
  template <typename V, typename CV> KOKKOS_INLINE_FUNCTION
  static void copy (V& d, const CV& s) {
    d[0] = s[0];
    d[1] = s[1];
    d[2] = s[2];
  }
  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void combine (const CV u, const CV v, const double a, V x) {
    const double oma = 1 - a;
    x[0] = oma*u[0] + a*v[0];
    x[1] = oma*u[1] + a*v[1];
    x[2] = oma*u[2] + a*v[2];
  }

  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void edge_normal (const CV a, const CV b, V en) {
    cross(a, b, en);
    normalize(en);
  }

  template <typename CV> KOKKOS_INLINE_FUNCTION
  static bool inside (const CV v, const CV a1, const CV a2, const CV n) {
    return dot_c_amb(n, v, a1) > 0 && dot_c_amb(n, v, a2) > 0;
  }

  /* Let
       n = edge normal
       c = edge point
       d = n' c
       v(a) = (1 - a) v1 + a v2.
     Solve n' v = d for a:
       a = (d - n' v1) / (n' (v2 - v1)).
     Then uvec(v(a)) is the intersection point on the unit sphere. Assume
     intersection exists. (Already filtered by 'inside'.)
  */
  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void intersect (const CV v1, const CV v2, const CV e1, const CV en,
                         V intersection) {
    double a; {
      const double
        num = dot_c_amb(en, e1, v1),
        den = dot_c_amb(en, v2, v1);
      a = num == 0 || den == 0 ? 0 : num/den;
      a = a < 0 ? 0 : a > 1 ? 1 : a;
    }
    // FP identity is sometimes useful, so let's enforce it.
    combine(v1, v2, a, intersection);
    normalize(intersection);
  }

  template <typename CV> KOKKOS_INLINE_FUNCTION
  static bool output (const CV v, int& no, Array2D<double>& vo) {
#ifdef SIKQ_DEBUG
    if (no >= vo.n()) {
      std::stringstream ss;
      ss << "output: No room in vo; vo.n() is " << vo.n() << " but no is "
         << no << "\n";
      error(ss.str().c_str());
    }
#endif
    if (no >= vo.n()) return false;
    vo(0,no) = v[0];
    vo(1,no) = v[1];
    vo(2,no) = v[2];
    ++no;
    return true;
  }

  //todo Handle non-convex case.
  // This uses a terrible formula, but it's just for testing.
  KOKKOS_INLINE_FUNCTION
  static double calc_area (const Array2D<const double>& v) {
    double area = 0;
    for (int i = 1; i < v.n() - 1; ++i) {
      const double a = calc_arc_length(v(0), v(i));
      const double b = calc_arc_length(v(i), v(i+1));
      const double c = calc_arc_length(v(i+1), v(0));
      const double s = 0.5*(a + b + c);
      const double d = (std::tan(0.5*s)*std::tan(0.5*(s-a))*
                        std::tan(0.5*(s-b))*std::tan(0.5*(s-c)));
      if (d <= 0) continue;
      area += 4*std::atan(std::sqrt(d));
    }
    return area;
  }
  template <typename CV> KOKKOS_INLINE_FUNCTION
  static double calc_arc_length (const CV a, const CV b) {
    const double d = dot(a, b);
    if (d >= 1) return 0;
    return acos(d);
  }
};

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
  Array2D<const double> p, nml;
  Array2D<const int> e, en;
};

// Generally not a user routine.
template <typename geo, typename CV> KOKKOS_INLINE_FUNCTION
bool clip_against_edge (
  // Input vertex list.
  const Array2D<const double>& vi, const int ni,
  // Output vertex list.
  Array2D<double>& vo, int& no,
  // The end points of the clip edge segment.
  const CV ce1, const CV ce2,
  // Clip edge's inward-facing normal.
  const CV cen)
{
  const double* s, * p;
  double intersection[3];
  no = 0;
  s = vi(ni-1);
  for (int j = 0; j < ni; ++j) {
    p = vi(j);
    if (geo::inside(p, ce1, ce2, cen)) {
      if (geo::inside(s, ce1, ce2, cen)) {
        if ( ! geo::output(p, no, vo)) return false;
      } else {
        geo::intersect(s, p, ce1, cen, intersection);
        if ( ! geo::output(intersection, no, vo)) return false;
        if ( ! geo::output(p, no, vo)) return false;
      }
    } else if (geo::inside(s, ce1, ce2, cen)) {
      geo::intersect(s, p, ce1, cen, intersection);
      if ( ! geo::output(intersection, no, vo)) return false;
    }
    s = p;
  }
  return true;
}

// Efficient user routine that uses the mesh data structure.
template <typename geo> KOKKOS_INLINE_FUNCTION
bool clip_against_poly (
  // Clip mesh. m.e(:,cp_e) is the element, and m.en(:,cp_e) is the
  // corresponding list of normal indices.
  const Mesh& m, const int cp_e,
  // A list of vertices describing the polygon to clip. The vertices must be in
  // a convention-determined order, such as CCW. vi(:,1:ni-1) are valid entries.
  const Array2D<const double>& vi, const int ni,
  // On output, vo(:,0:no-1) are vertices of the clipped polygon. no is 0 if
  // there is no intersection.
  Array2D<double>& vo, int& no,
  // Workspace. nvertwrk applies to both wrk and vo.n(). If nvertwrk is not
  // large enough, false is returned.
  double* const wrk, const int nvertwrk)
{
  Array2D<double> vo1(3, nvertwrk, wrk);
  int nos[] = { 0, 0 };
  Array2D<double>* vs[] = { &vo, &vo1 };

  const auto e = m.e(cp_e);
  const auto en = m.en(cp_e);

  auto nv = m.e.m(); // Number of vertices in clip polygon.
  while (e[nv-1] == -1) --nv;

  no = 0;
  if (nv % 2 == 0) {
    // Make sure the final vertex output list is in the caller's buffer.
    std::swap(vs[0], vs[1]);
    std::swap(nos[0], nos[1]);
  }

  if ( ! clip_against_edge<geo>(vi, ni, *vs[0], nos[0], m.p(e[0]), m.p(e[1]),
                                m.nml(en[0])))
    return false;
  if ( ! nos[0]) return true;

  for (int ie = 1, ielim = nv - 1; ; ++ie) {
    if ( ! clip_against_edge<geo>(*vs[0], nos[0], *vs[1], nos[1], m.p(e[ie]),
                                  m.p(e[(ie+1) % nv]), m.nml(en[ie])))
      return false;
    if ( ! nos[1]) return true;
    if (ie == ielim) break;
    std::swap(vs[0], vs[1]);
    std::swap(nos[0], nos[1]);
  }

  no = nos[1];
  return true;
}

// Not used for real stuff; just a convenient version for testing. In this
// version, clip_poly is a list of clip polygon vertices. This is instead of the
// mesh data structure.
template <typename geo> KOKKOS_INLINE_FUNCTION
bool clip_against_poly (
  // Clip polygon.
  const Array2D<const double>& clip_poly,
  // Clip polygon edges' inward-facing normals.
  const Array2D<const double>& clip_edge_normals,
  const Array2D<const double>& vi, const int ni,
  Array2D<double>& vo, int& no,
  double* const wrk, const int nvertwrk)
{
  Array2D<double> vo1(3, nvertwrk, wrk);
  int nos[] = { 0, 0 };
  Array2D<double>* vs[] = { &vo, &vo1 };

  no = 0;
  const auto nv = clip_poly.n();
  if (nv % 2 == 0) {
    // Make sure the final vertex output list is in the caller's buffer.
    std::swap(vs[0], vs[1]);
    std::swap(nos[0], nos[1]);
  }

  if ( ! clip_against_edge<geo>(vi, ni, *vs[0], nos[0], clip_poly(0), clip_poly(1),
                                clip_edge_normals(0)))
    return false;
  if ( ! nos[0]) return true;

  for (int ie = 1, ielim = nv - 1; ; ++ie) {
    if ( ! clip_against_edge<geo>(*vs[0], nos[0], *vs[1], nos[1], clip_poly(ie),
                                  clip_poly((ie+1) % nv), clip_edge_normals(ie)))
      return false;
    if ( ! nos[1]) return true;
    if (ie == ielim) break;
    std::swap(vs[0], vs[1]);
    std::swap(nos[0], nos[1]);
  }

  no = nos[1];
  return true;
}
} // namespace sh

// Octree for search.
class Octree {
public:
  typedef double BoundingBox[6];

  struct Options {
    // Do not go beyond max_depth depth, including the root and leaf. With this
    // constraint, try to go deep enough so that a leaf has no more than
    // max_nelem elements.
    int max_depth, max_nelem;
    Options () : max_depth(10), max_nelem(8) {}
  };

  // Bounding box for a cluster of points ps (possibly vertices).
  static void calc_bb (const Array2D<const double>& ps, const int np,
                       BoundingBox bb) {
    if (np == 0) return;
    for (int j = 0; j < 3; ++j)
      bb[j] = bb[j+3] = ps(j,0);
    for (int i = 1; i < np; ++i)
      for (int j = 0; j < 3; ++j) {
        bb[j] = std::min(bb[j], ps(j,i));
        bb[j+3] = std::max(bb[j+3], ps(j,i));
      }
  }

  static void calc_bb (const Array2D<const double>& ps, BoundingBox bb) {
    calc_bb(ps, ps.n(), bb);
  }

  template <typename CIV, typename V>
  static void calc_bb (const Array2D<const double>& p, const CIV& e,
                       const int ne, V ebb) {
    for (int j = 0; j < 3; ++j)
      ebb[j] = ebb[j+3] = p(j, e[0]);
    for (int i = 1; i < ne; ++i) {
      if (e[i] == -1) break;
      for (int j = 0; j < 3; ++j) {
        ebb[j] = ko::min(ebb[j], p(j, e[i]));
        ebb[j+3] = ko::max(ebb[j+3], p(j, e[i]));
      }
    }
  }

  static void calc_bb (const Array2D<const double>& p, const Array2D<const int>& e,
                       Array2D<double>& ebbs) {
    assert(ebbs.n() == e.n());
    for (int k = 0; k < e.n(); ++k)
      calc_bb(p, e(k), e.m(), ebbs(k));
  }

  // p is a 3xNp array of points. e is a KxNe array of elements. An entry <0 is
  // ignored. All <0 entries must be at the end of an element's list.
  Octree (const Array2D<const double>& p, const Array2D<const int>& e,
          const Options& o) {
    init(p, e, o);
  }
  Octree (const Array2D<const double>& p, const Array2D<const int>& e) {
    Options o;
    init(p, e, o);
  }

  // Apply f to every element in leaf nodes with which bb overlaps. f must have
  // function
  //     void operator(const int element_index).
  // element_index indexes e.
  template <typename CV, typename Functor> KOKKOS_INLINE_FUNCTION
  void apply (const CV bb, Functor& f) const {
    if (nodes_.n() == 0) {
      for (int i = 0; i < offset_[1]; ++i)
        f(elems_[i]);
      return;
    }
    apply_r(0, bb_, bb, f);
  }

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

  class IntList {
    int* const buf_;
    int i_;
  public:
    IntList (int* const buf) : buf_(buf), i_(0) {}
    void reset () { i_ = 0; }
    void push (const int& i) { buf_[i_++] = i; }
    int* data () { return buf_; }
    int n () const { return i_; }
    const int& operator[] (const int& i) const { return buf_[i]; }
  };

  class DynIntList {
    std::vector<int> buf_;
  public:
    DynIntList () {}
    void push (const int& i) { buf_.push_back(i); }
    int& back () { return buf_.back(); }
    int& operator[] (const size_t i) {
      if (i >= buf_.size())
        buf_.resize(i+1);
      return buf_[i];
    }
    int n () const { return static_cast<int>(buf_.size()); }
    const int* data () const { return buf_.data(); }
  };

  class Nodes {
    std::vector<int> buf_;
  public:
    int n () const { return static_cast<int>(buf_.size()) >> 3; }
    const int* data () const { return buf_.data(); }
    int& operator() (const int& r, const int& c) {
      const size_t ec = (c+1) << 3;
      if (ec >= buf_.size())
        buf_.resize(ec);
      assert(((c << 3) + r) >= 0);
      assert(((c << 3) + r) < (int) buf_.size());
      return buf_[(c << 3) + r];
    }
  };

  void init (const Array2D<const double>& p, const Array2D<const int>& e,
             const Options& o) {
    if (e.n() == 0) return;
    // Get OT's bounding box.
    calc_bb(p, bb_);
    // Get elements' bounding boxes.
    Array2D<double> ebbs(6, e.n());
    calc_bb(p, e, ebbs);
    // Static element lists for work. Each level has active work space.
    std::vector<int> buf((o.max_depth - 1)*e.n());
    IntList es(buf.data()), wrk(buf.data() + e.n());
    for (int i = 0; i < e.n(); ++i)
      es.push(i);
    // Dynamic element lists.
    DynIntList offset, elems;
    offset[0] = 0;
    // Dynamic node data structure.
    Nodes nodes;
    // Recurse. We don't care about the return value. If it's 0 and nodes.n() ==
    // 0, we'll detect as much in 'apply'.
    init_r(1, bb_, ebbs, o, es, wrk, offset, elems, nodes);
    // Build the static data structures.
    if (elems.n() == 0) return;
    offset_.reset(offset.n());
    elems_.reset(elems.n());    
    memcpy(offset_.data(), offset.data(), offset.n() * sizeof(*offset_.data()));
    memcpy(elems_.data(), elems.data(), elems.n() * sizeof(*offset_.data()));
    if (nodes.n() == 0) return;
    nodes_.reset(8, nodes.n());
    memcpy(nodes_.data(), nodes.data(), (nodes.n() << 3) * sizeof(*offset_.data()));
    // Move them to the device.
    nodes_.modify(); nodes_.device().sync();
    offset_.modify(); offset_.device().sync();
    elems_.modify(); elems_.device().sync();
  }

  int init_r (const int depth, // Tree's depth at this point, including root.
              const BoundingBox& nbb, // My bounding box.
              const Array2D<const double>& ebbs, // All elements' bounding boxes.
              const Options& o, // Options controlling construct of the tree.
              IntList& es, // List of elements in my bounding box.
              IntList& wrk, // Work space to store working element lists.
              DynIntList& offset, // Offsets into elems.
              DynIntList& elems, // Elements belonging to leaf nodes.
              Nodes& nodes) // Dynamic nodes data structure.
  {
    const int my_idx = nodes.n(); // My node index.
    // Decide what to do.
    if (es.n() == 0) {
      // I have no elements, so return 0 to indicate I'm a leaf node containing
      // nothing.
      return 0;
    } else if (es.n() <= o.max_nelem || depth == o.max_depth) {
      // I'm a leaf node with elements. Store my list of elements and return the
      // storage location.
      const int os = offset.back();
      offset.push(os + es.n());
      for (int i = 0, n = es.n(); i < n; ++i)
        elems[os + i] = es[i];
      return 1 - offset.n();
    } else {
      // I'm not a leaf node.
      nodes(0, my_idx) = 0; // Insert myself into the nodes array.
      for (int ic = 0; ic < 8; ++ic) {
        BoundingBox child_bb;
        fill_child_bb(nbb, ic, child_bb);
        // Find the elements that are in this child's bb.
        IntList ces(wrk.data());
        for (int i = 0, n = es.n(); i < n; ++i)
          if (do_bb_overlap(child_bb, ebbs(es[i])))
            ces.push(es[i]);
        // Create some work space.
        IntList cwrk(wrk.data() + ces.n());
        // Recurse.
        const int child_idx = init_r(depth+1, child_bb, ebbs, o, ces, cwrk,
                                     offset, elems, nodes);
        nodes(ic, my_idx) = child_idx;
      }
      return my_idx;
    }
  }

  // Using parent bb p, fill child bb c, with child_idx in 0:7.
  KOKKOS_INLINE_FUNCTION
  static void fill_child_bb (const BoundingBox& p, const int& child_idx,
                             BoundingBox& c) {
    const double m[] = { 0.5*(p[0] + p[3]),
                         0.5*(p[1] + p[4]),
                         0.5*(p[2] + p[5]) };
    switch (child_idx) {
    case 0: c[0] = p[0]; c[1] = p[1]; c[2] = p[2]; c[3] = m[0]; c[4] = m[1]; c[5] = m[2]; break;
    case 1: c[0] = m[0]; c[1] = p[1]; c[2] = p[2]; c[3] = p[3]; c[4] = m[1]; c[5] = m[2]; break;
    case 2: c[0] = m[0]; c[1] = m[1]; c[2] = p[2]; c[3] = p[3]; c[4] = p[4]; c[5] = m[2]; break;
    case 3: c[0] = p[0]; c[1] = m[1]; c[2] = p[2]; c[3] = m[0]; c[4] = p[4]; c[5] = m[2]; break;
    case 4: c[0] = p[0]; c[1] = p[1]; c[2] = m[2]; c[3] = m[0]; c[4] = m[1]; c[5] = p[5]; break;
    case 5: c[0] = m[0]; c[1] = p[1]; c[2] = m[2]; c[3] = p[3]; c[4] = m[1]; c[5] = p[5]; break;
    case 6: c[0] = m[0]; c[1] = m[1]; c[2] = m[2]; c[3] = p[3]; c[4] = p[4]; c[5] = p[5]; break;
    case 7: c[0] = p[0]; c[1] = m[1]; c[2] = m[2]; c[3] = m[0]; c[4] = p[4]; c[5] = p[5]; break;
    default:
      // impossible
      error("fill_child_bb: The impossible has happened.");
    }
  }

  // Do bounding boxes a and b overlap?
  KOKKOS_INLINE_FUNCTION
  static bool do_bb_overlap (const BoundingBox a, const BoundingBox b) {
    for (int i = 0; i < 3; ++i)
      if ( ! do_lines_overlap(a[i], a[i+3], b[i], b[i+3]))
        return false;
    return true;
  }

  KOKKOS_INLINE_FUNCTION
  static bool do_lines_overlap (const double& a1, const double& a2,
                                const double& b1, const double& b2) {
    return ! (a2 < b1 || a1 > b2);
  }

  template <typename CV, typename Functor> KOKKOS_INLINE_FUNCTION
  void apply_r (const int ni, const BoundingBox& nbb, const CV bb,
                Functor& f) const {
    for (int i = 0; i < 8; ++i) {
      BoundingBox child_bb;
      fill_child_bb(nbb, i, child_bb);
      if ( ! do_bb_overlap(child_bb, bb)) continue;
      int e = nodes_(i,ni);
      if (e > 0)
        apply_r(e, child_bb, bb, f);
      else if (e < 0) {
        e = std::abs(e + 1);
        for (int k = offset_[e]; k < offset_[e+1]; ++k)
          f(elems_[k]);
      }
    }
  }
};

namespace test {
static constexpr int max_nvert = 20;

// In practice, we want to form high-quality normals using information about the
// mesh, such as that it is a CS mesh. For testing, form the normals from edge
// vertices. (This leads to increasing cancellation error with mesh refinement.)
template <typename geo>
void fill_normals (sh::Mesh& m) {
  // Count number of edges.
  int ne = 0;
  for (auto ip = zero(m.e); ip < m.e.n(); ++ip)
    for (auto iv = zero(m.e); iv < m.e.m(); ++iv)
      if (m.e(iv,ip) == -1) break; else ++ne;
  // Fill.
  Array2D<int> en(m.e.m(), m.e.n());
  en.set(-1);
  Array2D<double> nml(3, ne);
  int ie = 0;
  for (auto ip = zero(m.e); ip < m.e.n(); ++ip)
    for (auto iv = zero(m.e); iv < m.e.m(); ++iv)
      if (m.e(iv,ip) == -1)
        break;
      else {
        // Somewhat complicated next node index.
        const int iv_next = (iv+1 == m.e.m() ? 0 :
                             (m.e(iv+1,ip) == -1 ? 0 : iv+1));
        geo::edge_normal(m.p(m.e(iv, ip)), m.p(m.e(iv_next, ip)), nml(ie));
        en(iv,ip) = ie;
        ++ie;
      }
  m.en = en;
  m.nml = nml;
}

// Used in Octree::apply to gather a set of possibly intersecting polygons.
struct OTSearchFunctor {
  std::set<int> hits;
  KOKKOS_INLINE_FUNCTION void operator() (const int i) { hits.insert(i); }
};

// Find the area of the overlapping part of two meshes by summing over the areas
// of the common refinement polygons. Obviously a silly thing to do, but a good
// test and demonstration problem.
template <typename geo>
class TestAreaOTFunctor {
  // Mesh against which to clip. ("Eulerian mesh".)
  sh::Mesh cm;
  // Mesh of clipped polygons. ("Departure mesh".)
  const Array2D<const double> p; // 3 x #verts array of polygon vertices.
  const Array2D<const int> e;    // Array of polygons. e(:,k) is the k'th polygon.
  // Already initialized octree used to search for possibly intersecting
  // polygons.
  Octree ot;

public:
  typedef double value_type;

  TestAreaOTFunctor (const sh::Mesh& cm, const Array2D<const double>& p,
                     const Array2D<const int>& e, const Octree& ot)
    : cm(cm), p(p), e(e), ot(ot)
  {}
  
  // k indexes (p,e).
  KOKKOS_INLINE_FUNCTION void operator() (const int k, double& area) const {
    // Clipped element bounding box.
    double ebb[6];
    Octree::calc_bb(p, e(k), e.m(), ebb);
    // Get list of possible overlaps.
    OTSearchFunctor f;
    ot.apply(ebb, f);
    // In and out vertex lists.
    double buf[6*max_nvert];
    Array2D<double>
      vi(3, max_nvert, buf),
      vo(3, max_nvert, buf + 3*max_nvert);
    int ni, no;
    // Workspace.
    double wrk[3*max_nvert];
    // Area of all overlapping regions.
    double a = 0;
    for (const auto icp : f.hits) {
      // Create the polygon to be clipped.
      ni = 0;
      for (int i = 0; i < e.m(); ++i) {
        if (e(i,k) == -1) break;
        copy(vi(i), p(e(i,k)), 3);
        ++ni;
      }
      sh::clip_against_poly<geo>(cm, icp, vi, ni, vo, no, wrk, max_nvert);
      if (no) {
        // A non-0 intersection was found. Accumulate the area.
        a += geo::calc_area(Array2D<const double>(vo.m(), no, vo.data()));
      }
    }
    // Add our area to the reduction.
    area += a;
  }
};

template <typename geo>
double test_area_ot (const Array2D<const double>& cp, const Array2D<const int>& ce,
                     const Array2D<const double>& p, const Array2D<const int>& e) {
  // Clip mesh and edge normal calculation. (In practice, we'd like to use
  // higher-quality edge normals.)
  sh::Mesh cm; cm.p = cp; cm.e = ce;
  fill_normals<geo>(cm);

  double et[2];
  auto t = tic();
  // Build an octree over the clip mesh.
  Octree ot(cp, ce);
  et[0] = toc(t);

  // Compute the area in a silly way to test search and interesection.
  t = tic();
  double area = 0;
  ko::parallel_reduce(e.n(), TestAreaOTFunctor<geo>(cm, p, e, ot), area);
  et[1] = toc(t);
  print_times("test_area_ot", et, 2);
  return area;
}
} // namespace test
} // namespace siqk

#endif // INCLUDE_SIK_HPP
