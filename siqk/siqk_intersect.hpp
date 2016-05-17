#ifndef INCLUDE_SIQK_INTERSECT_HPP
#define INCLUDE_SIQK_INTERSECT_HPP

#include "siqk.hpp"
#include "siqk_quadrature.hpp"

namespace siqk {
struct PlaneGeometry {
  template <typename V> KOKKOS_INLINE_FUNCTION
  static void scale (const Real& a, V v) {
    v[0] *= a; v[1] *= a;
  }
  template <typename CV> KOKKOS_INLINE_FUNCTION
  static Real dot_c_amb (const CV c, const CV a, const CV b) {
    return c[0]*(a[0] - b[0]) + c[1]*(a[1] - b[1]);
  }
  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void combine (const CV u, const CV v, const Real& a, V x) {
    const Real& oma = 1 - a;
    x[0] = oma*u[0] + a*v[0];
    x[1] = oma*u[1] + a*v[1];
  }

  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void edge_normal (const CV e1, const CV e2, V en) {
    en[0] = e1[1] - e2[1];
    en[1] = e2[0] - e1[0];
  }

  template <typename CV> KOKKOS_INLINE_FUNCTION
  static bool inside (const CV v, const CV e1, const CV en) {
    return dot_c_amb(en, v, e1) >= 0;
  }

  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void intersect (const CV v1, const CV v2, const CV e1, const CV en,
                         V intersection) {
    const Real& a = dot_c_amb(en, e1, v1) / dot_c_amb(en, v2, v1);
    combine(v1, v2, a, intersection);
  }

  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static bool output (const CV v, Int& no, const V& vo) {
#ifdef SIQK_DEBUG
    if (no >= nslices(vo)) {
      std::stringstream ss;
      ss << "output: No room in vo; vo.n() is " << vo.n() << " but no is "
         << no << "\n";
      error(ss.str().c_str());
    }
#endif
    if (no >= nslices(vo)) return false;
    vo(no,0) = v[0];
    vo(no,1) = v[1];
    ++no;
    return true;
  }

  //todo Handle non-convex case.
  template <typename CV3s>
  KOKKOS_INLINE_FUNCTION
  static Real calc_area (const CV3s& v, const Int n) {
    Real area = 0;
    for (Int i = 1, ilim = n - 1; i < ilim; ++i) {
      Real v1[2], v2[2];
      v1[0] = v(i,0) - v(0,0);
      v1[1] = v(i,1) - v(0,1);
      v2[0] = v(i+1,0) - v(0,0);
      v2[1] = v(i+1,1) - v(0,1);
      const Real a = v1[0]*v2[1] - v1[1]*v2[0];
      area += a;
    }
    return 0.5*area;
  }

  template <typename CV3s>
  KOKKOS_INLINE_FUNCTION
  static Real calc_area_formula (const CV3s& v, const Int n) {
    return calc_area(v, n);
  }
};

// All inputs and outputs are relative to the unit-radius sphere.
struct SphereGeometry {
  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void cross (const CV a, const CV b, V c) {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
  }
  template <typename CV> KOKKOS_INLINE_FUNCTION
  static Real dot (const CV a, const CV b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  }
  template <typename CV> KOKKOS_INLINE_FUNCTION
  static Real norm2 (const CV v) {
    return dot(v, v);
  }
  template <typename V> KOKKOS_INLINE_FUNCTION
  static void scale (const Real& a, V v) {
    v[0] *= a; v[1] *= a; v[2] *= a;
  }
  template <typename V> KOKKOS_INLINE_FUNCTION
  static void normalize (V v) {
    scale(1.0/std::sqrt(norm2(v)), v);
  }
  template <typename CV> KOKKOS_INLINE_FUNCTION
  static Real dot_c_amb (const CV c, const CV a, const CV b) {
    return c[0]*(a[0] - b[0]) + c[1]*(a[1] - b[1]) + c[2]*(a[2] - b[2]);
  }
  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void axpy (const Real& a, const CV& x, V& y) {
    y[0] += a*x[0];
    y[1] += a*x[1];
    y[2] += a*x[2];
  }
  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void axpbyz (const Real& a, const CV& x, const Real& b, const CV& y,
                      V& z) {
    z[0] = a*x[0] + b*y[0];
    z[1] = a*x[1] + b*y[1];
    z[2] = a*x[2] + b*y[2];
  }
  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void combine (const CV& u, const CV& v, const Real& a, V& x) {
    const Real& oma = 1 - a;
    x[0] = oma*u[0] + a*v[0];
    x[1] = oma*u[1] + a*v[1];
    x[2] = oma*u[2] + a*v[2];
  }

  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void edge_normal (const CV a, const CV b, V en) {
    cross(a, b, en);
    normalize(en);
  }

  // Is v inside the line anchored at a with inward-facing normal n?
  template <typename CV> KOKKOS_INLINE_FUNCTION
  static bool inside (const CV& v, const CV& a, const CV& n) {
    return dot_c_amb(n, v, a) >= 0;
  }

  /* Let
       en = edge normal
       e1 = edge starting point
       d = en' e1
       v(a) = (1 - a) v1 + a v2.
     Solve n' v = d for a:
       a = (en' (e1 - v1)) / (en' (v2 - v1)).
     Then uvec(v(a)) is the intersection point on the unit sphere. Assume
     intersection exists. (Already filtered by 'inside'.)
  */
  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void intersect (const CV v1, const CV v2, const CV e1, const CV en,
                         V intersection) {
    const Real a = dot_c_amb(en, e1, v1) / dot_c_amb(en, v2, v1);
    combine(v1, v2, a, intersection);
    normalize(intersection);
  }

  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static bool output (const CV v, Int& no, V& vo) {
    if (no >= nslices(vo)) return false;
    vo(no,0) = v[0];
    vo(no,1) = v[1];
    vo(no,2) = v[2];
    ++no;
    return true;
  }

  //todo Handle non-convex case.
  // This uses a terrible formula, but it's just for testing.
  template <typename CV3s>
  KOKKOS_INLINE_FUNCTION
  static Real calc_area_formula (const CV3s& v, const Int n) {
    Real area = 0;
    for (Int i = 1, ilim = n - 1; i < ilim; ++i) {
      const Real a = calc_arc_length(slice(v,0), slice(v,i));
      const Real b = calc_arc_length(slice(v,i), slice(v,i+1));
      const Real c = calc_arc_length(slice(v,i+1), slice(v,0));
      const Real s = 0.5*(a + b + c);
      const Real d = (std::tan(0.5*s)*std::tan(0.5*(s-a))*
                      std::tan(0.5*(s-b))*std::tan(0.5*(s-c)));
      if (d <= 0) continue;
      area += 4*std::atan(std::sqrt(d));
    }
    return area;
  }
  template <typename CV> KOKKOS_INLINE_FUNCTION
  static Real calc_arc_length (const CV a, const CV b) {
    const Real d = dot(a, b);
    if (d >= 1) return 0;
    return acos(d);
  }

  template <typename CV3s, typename ExeSpace = ko::DefaultExecutionSpace>
  KOKKOS_INLINE_FUNCTION
  static Real calc_area (const CV3s& v, const Int n) {
    Real area = 0;
    for (Int i = 1, ilim = n - 1; i < ilim; ++i) {
      Real a = 0;
      RawConstVec3s coord;
      RawConstArray weight;
      quadrature::get_coef<ExeSpace>(4, coord, weight);
      for (Int k = 0, klim = nslices(coord); k < klim; ++k) {
        const Real jac = calc_tri_jacobian(slice(v,0), slice(v,i), slice(v,i+1),
                                           slice(coord, k));
        a += weight[k]*jac;
      }
      area += 0.5*a;
    }
    return area;
  }
  template <typename CV, typename CA>
  KOKKOS_INLINE_FUNCTION
  static Real calc_tri_jacobian (const CV& v1, const CV& v2, const CV& v3,
                                 const CA& alpha) {
    // V(:,i) is vertex i of the spherical triangle on the unit sphere. The
    // coefs
    //     alpha = [a1, a2, 1 - a1 - a2]'
    //           = [1 0; 0 1; -1 -1] [a1, a2]'
    //           = alpha_a a
    // (barycentric coords) give the location
    //     v = V alpha
    // on the planar triangle, and u = uvec(v) is the point on the unit sphere.
    //   For a planar tri in 3D, the jacobian is
    //     v_a = v_alpha alpha_a
    //         = V [1 0; 0 1; -1 -1]
    //     J = norm(cross(v_a(:,1), v_a(:,2))).
    // For a spherical tri with the same vertices,
    //     u = v/(v' v)
    //     u_a = u_alpha alpha_a
    //         = (v'v)^{-1/2} (I - u u') V alpha_a
    //         = (v'v)^{-1/2} (I - u u') v_a
    //     J = norm(cross(u_a(:,1), u_a(:,2))).
    Real u[3] = {0};
    axpy(alpha[0], v1, u);
    axpy(alpha[1], v2, u);
    axpy(alpha[2], v3, u);
    const auto oovn = 1/std::sqrt(norm2(u));
    scale(oovn, u);
    Real u_a[2][3];
    axpbyz(1, v1, -1, v3, u_a[0]);
    axpbyz(1, v2, -1, v3, u_a[1]);
    for (int i = 0; i < 2; ++i) {
      axpy(-dot(u, u_a[i]), u, u_a[i]);
      scale(oovn, u_a[i]);
    }
    cross(u_a[0], u_a[1], u);
    return std::sqrt(norm2(u));
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
template <typename ES = ko::DefaultExecutionSpace>
struct Mesh {
  typename InExeSpace<ConstVec3s, ES>::type p, nml;
  typename InExeSpace<ConstIdxs, ES>::type e, en;

  Mesh () {}

  Mesh (const Mesh<ko::HostSpace>& m) {
    typename InExeSpace<Vec3s, ES>::type tp, tnml;
    typename InExeSpace<Idxs, ES>::type te, ten;
    resize_and_copy(tp, m.p); p = tp;
    resize_and_copy(tnml, m.nml); nml = tnml;
    resize_and_copy(te, m.e); e = te;
    resize_and_copy(ten, m.en); en = ten;
  }
};

// Generally not a user routine.
template <typename geo, typename CV3s, typename V3s, typename CV>
KOKKOS_INLINE_FUNCTION
bool clip_against_edge (
  // Input vertex list.
  const CV3s& vi, const Int ni,
  // Output vertex list.
  V3s& vo, Int& no,
  // One point of the clip edge.
  const CV ce1,
  // Clip edge's inward-facing normal.
  const CV cen)
{
  Real intersection[3];
  no = 0;
  auto s = const_slice(vi, ni-1);
  for (Int j = 0; j < ni; ++j) {
    auto p = const_slice(vi,j);
    if (geo::inside(p, ce1, cen)) {
      if (geo::inside(s, ce1, cen)) {
        if ( ! geo::output(p, no, vo)) return false;
      } else {
        geo::intersect(s, p, ce1, cen, intersection);
        if ( ! geo::output(intersection, no, vo)) return false;
        if ( ! geo::output(p, no, vo)) return false;
      }
    } else if (geo::inside(s, ce1, cen)) {
      geo::intersect(s, p, ce1, cen, intersection);
      if ( ! geo::output(intersection, no, vo)) return false;
    }
    s = p;
  }
  return true;
}

// Efficient user routine that uses the mesh data structure.
//todo An optimization would be to have 2 clip_against_edge routines. One would
// handle the special case of the first vertex list being in (p,e) format.
template <typename geo, typename MeshT, typename CV3s, typename V3s>
KOKKOS_INLINE_FUNCTION
bool clip_against_poly (
  // Clip mesh. m.e(:,cp_e) is the element, and m.en(:,cp_e) is the
  // corresponding list of normal indices.
  const MeshT& m, const Int cp_e,
  // A list of vertices describing the polygon to clip. The vertices must be in
  // a convention-determined order, such as CCW. vi(:,1:ni-1) are valid entries.
  const CV3s& vi, const Int ni,
  // On output, vo(:,0:no-1) are vertices of the clipped polygon. no is 0 if
  // there is no intersection.
  V3s& vo, Int& no,
  // Workspace. Both vo and wrk must be large enough to hold all generated
  // vertices. If they are not, false is returned.
  V3s& wrk)
{
  Int nos[] = { 0, 0 };
  V3s* vs[] = { &vo, &wrk };

  const auto e = slice(m.e, cp_e);
  const auto en = slice(m.en, cp_e);

  auto nv = szslice(m.e); // Number of vertices in clip polygon.
  while (e[nv-1] == -1) --nv;

  no = 0;
  if (nv % 2 == 0) {
    // Make sure the final vertex output list is in the caller's buffer.
    swap(vs[0], vs[1]);
    swap(nos[0], nos[1]);
  }

  if ( ! clip_against_edge<geo>(vi, ni, *vs[0], nos[0], const_slice(m.p, e[0]),
                                const_slice(m.nml, en[0])))
    return false;
  if ( ! nos[0]) return true;

  for (Int ie = 1, ielim = nv - 1; ; ++ie) {
    if ( ! clip_against_edge<geo>(*vs[0], nos[0], *vs[1], nos[1],
                                  const_slice(m.p, e[ie]),
                                  const_slice(m.nml, en[ie])))
      return false;
    if ( ! nos[1]) return true;
    if (ie == ielim) break;
    swap(vs[0], vs[1]);
    swap(nos[0], nos[1]);
  }

  no = nos[1];
  return true;
}

// Not used for real stuff; just a convenient version for testing. In this
// version, clip_poly is a list of clip polygon vertices. This is instead of the
// mesh data structure.
template <typename geo, typename CV3s, typename V3s>
KOKKOS_INLINE_FUNCTION
bool clip_against_poly (
  // Clip polygon.
  const CV3s& clip_poly,
  // Clip polygon edges' inward-facing normals.
  const CV3s& clip_edge_normals,
  const CV3s& vi, const Int ni,
  V3s& vo, Int& no,
  V3s& wrk)
{
  Int nos[] = { 0, 0 };
  V3s* vs[] = { &vo, &wrk };

  no = 0;
  if (nslices(clip_poly) % 2 == 0) {
    // Make sure the final vertex output list is in the caller's buffer.
    swap(vs[0], vs[1]);
    swap(nos[0], nos[1]);
  }

  if ( ! clip_against_edge<geo>(vi, ni, *vs[0], nos[0],
                                const_slice(clip_poly, 0),
                                const_slice(clip_edge_normals, 0)))
    return false;
  if ( ! nos[0]) return true;

  for (Int ie = 1, ielim = nslices(clip_poly) - 1; ; ++ie) {
    if ( ! clip_against_edge<geo>(*vs[0], nos[0], *vs[1], nos[1],
                                  const_slice(clip_poly, ie),
                                  const_slice(clip_edge_normals, ie)))
      return false;
    if ( ! nos[1]) return true;
    if (ie == ielim) break;
    swap(vs[0], vs[1]);
    swap(nos[0], nos[1]);
  }

  no = nos[1];
  return true;
}
} // namespace sh

// Oct-tree. Might do something else better suited to the sphere later.
template <Int max_depth_ = 10>
class Octree {
public:
  enum { max_depth = max_depth_ };
  typedef Real BoundingBox[6];

  struct Options {
    // Do not go beyond max_depth_ depth, including the root and leaf. With this
    // constraInt, try to go deep enough so that a leaf has no more than
    // max_nelem elements.
    Int max_nelem;
    Options () : max_nelem(8) {}
  };

  // Bounding box for a cluster of points ps (possibly vertices).
  //todo kernelize
  template <typename CV3s>
  static void calc_bb (const CV3s& ps, const Int np, BoundingBox bb) {
    if (np == 0) return;
    for (Int j = 0; j < 3; ++j)
      bb[j] = bb[j+3] = ps(0,j);
    for (Int i = 1; i < np; ++i)
      for (Int j = 0; j < 3; ++j) {
        bb[j] = min(bb[j], ps(i,j));
        bb[j+3] = max(bb[j+3], ps(i,j));
      }
  }

  template <typename CV3s>
  static void calc_bb (const CV3s& ps, BoundingBox bb) {
    calc_bb(ps, nslices(ps), bb);
  }

  template <typename CV3s, typename CIV, typename V>
  KOKKOS_INLINE_FUNCTION
  static void calc_bb (const CV3s& p, const CIV& e, const Int ne, V ebb) {
    for (Int j = 0; j < 3; ++j)
      ebb[j] = ebb[j+3] = p(e[0], j);
    for (Int i = 1; i < ne; ++i) {
      if (e[i] == -1) break;
      for (Int j = 0; j < 3; ++j) {
        ebb[j] = min(ebb[j], p(e[i], j));
        ebb[j+3] = max(ebb[j+3], p(e[i], j));
      }
    }
  }

  //todo kernelize
  template <typename CV3s, typename CIs, typename V6s>
  static void calc_bb (const CV3s& p, const CIs& e, V6s& ebbs) {
    assert(nslices(ebbs) == nslices(e));
    for (Int k = 0, klim = nslices(e); k < klim; ++k)
      calc_bb(p, slice(e, k), szslice(e), slice(ebbs, k));
  }

  // p is a 3xNp array of points. e is a KxNe array of elements. An entry <0 is
  // ignored. All <0 entries must be at the end of an element's list.
  Octree (const ConstVec3s::HostMirror& p, const ConstIdxs::HostMirror& e,
          const Options& o) {
    init(p, e, o);
  }
  Octree (const ConstVec3s::HostMirror& p, const ConstIdxs::HostMirror& e) {
    Options o;
    init(p, e, o);
  }

  // Apply f to every element in leaf nodes with which bb overlaps. f must have
  // function
  //     void operator(const Int element).
  template <typename CV, typename Functor>
  KOKKOS_INLINE_FUNCTION
  void apply (const CV bb, Functor& f) const {
    if (nslices(nodes_) == 0) {
      for (Int i = 0; i < offsets_[1]; ++i)
        f(elems_[i]);
      return;
    }
#ifdef SIQK_NONRECURSIVE
    // Non-recursive impl.
    {
      // Stack.
      Real snbb[8*max_depth_];
      Int sni[max_depth_], si[max_depth_];
      Int sp = 0;
      // Args for top-level call.
      copy(snbb, bb_, 8);
      sni[sp] = 0;
      si[sp] = 0;
      while (sp >= 0) {
        // Get stack frame's (nbb, ni, current i) values.
        const Int i = si[sp];
        if (i == 8) {
          --sp;
          continue;
        }
        // Increment stored value of i for next iteration. Current value is
        // stored in 'i' above.
        ++si[sp];
        const Int ni = sni[sp];
        const Real* const nbb = snbb + 8*sp;
        // Can use the next stack frame's bb space for a child bb.
        Real* const child_bb = snbb + 8*(sp+1);
        fill_child_bb(nbb, i, child_bb);
        if ( ! do_bb_overlap(child_bb, bb)) continue;
        Int e = nodes_(ni,i);
        if (e < 0) {
          // Leaf, so apply functor to each element.
          e = std::abs(e + 1);
          for (Int k = offsets_[e]; k < offsets_[e+1]; ++k)
            f(elems_[k]);
        } else if (e > 0) {
          // Recurse.
          ++sp;
          sni[sp] = e;
          si[sp] = 0;
        }
      }
    }
#else
    apply_r(0, bb_, bb, f);
#endif
  }

private:
  /* Each node in the oct-tree contains 8 integers, stored in 'nodes'.

     >0 is an index Into 'nodes', pointing to a child node.

     A <=0 entry in 'nodes' indicates a leaf node. If 0, there are no elements
     in the leaf. If <0, the negative of the entry minus 1 is the index of an
     offset array indexing 'elems'.

     Each segment of 'elems' contains a list of element indices covered by a
     leaf node. Element indices refer to the list of elements the caller
     provides during oct-tree construction.
  */

  // Static data structures holding the completed octree.
  //   nodes(:,i) is a list. The list includes children of node i (>0) and leaf
  // node data (<=0).
  //todo Make these const once ready to do full GPU stuff.
  Nodes nodes_;
  // A leaf node corresponding to -k covers elements
  //     elems[offset[k] : offset[k]-1].
  ko::View<int*> offsets_, elems_;
  // Root node's bounding box.
  BoundingBox bb_;

  // Dynamic data structures for construction phase.
  class IntList {
    Int* const buf_;
    Int i_;
  public:
    IntList (Int* const buf) : buf_(buf), i_(0) {}
    void reset () { i_ = 0; }
    void push (const Int& i) { buf_[i_++] = i; }
    Int* data () { return buf_; }
    Int n () const { return i_; }
    const Int& operator[] (const Int& i) const { return buf_[i]; }
  };

  class DynIntList {
    std::vector<Int> buf_;
  public:
    DynIntList () {}
    void push (const Int& i) { buf_.push_back(i); }
    Int& back () { return buf_.back(); }
    Int& operator[] (const size_t i) {
      if (i >= buf_.size())
        buf_.resize(i+1);
      return buf_[i];
    }
    const Int& operator[] (const size_t i) const { return buf_[i]; }
    Int n () const { return static_cast<Int>(buf_.size()); }
    const Int* data () const { return buf_.data(); }
  };

  // Opposite index slot convention.
  class DynNodes {
    std::vector<Int> buf_;
  public:
    Int n () const { return static_cast<Int>(buf_.size()) >> 3; }
    const Int* data () const { return buf_.data(); }
    Int& operator() (const Int& r, const Int& c) {
      const size_t ec = (c+1) << 3;
      if (ec >= buf_.size())
        buf_.resize(ec);
      return const_cast<Int&>(
        const_cast<const DynNodes*>(this)->operator()(r, c));
    }
    const Int& operator() (const Int& r, const Int& c) const {
      assert(((c << 3) + r) >= 0);
      assert(((c << 3) + r) < (Int) buf_.size());
      return buf_[(c << 3) + r];
    }
  };

  void init (const ConstVec3s::HostMirror& p, const ConstIdxs::HostMirror& e,
             const Options& o) {
    if (nslices(e) == 0) return;
    // Get OT's bounding box.
    calc_bb(p, bb_);
    // Get elements' bounding boxes.
    Vec6s::HostMirror ebbs("ebbs", nslices(e), 6);
    calc_bb(p, e, ebbs);
    // Static element lists for work. Each level has active work space.
    std::vector<Int> buf(max_depth_*nslices(e));
    IntList es(buf.data()), wrk(buf.data() + nslices(e));
    for (Int i = 0, ilim = nslices(e); i < ilim; ++i)
      es.push(i);
    // Dynamic element lists.
    DynIntList offsets, elems;
    offsets[0] = 0;
    // Dynamic node data structure.
    DynNodes nodes;
    // Recurse. We don't care about the return value. If it's 0 and nodes.n() ==
    // 0, we'll detect as much in 'apply'.
    init_r(1, bb_, ebbs, o, es, wrk, offsets, elems, nodes);
    // Build the static data structures.
    if (elems.n() == 0) return;
    init_static_ds(nodes, offsets, elems);
  }

  Int init_r (const Int depth, // Tree's depth at this point, including root.
              const BoundingBox& nbb, // My bounding box.
              const ConstVec6s::HostMirror& ebbs, // All elements' bounding boxes.
              const Options& o, // Options controlling construct of the tree.
              IntList& es, // List of elements in my bounding box.
              IntList& wrk, // Work space to store working element lists.
              DynIntList& offsets, // Offsetss Into elems.
              DynIntList& elems, // Elements belonging to leaf nodes.
              DynNodes& nodes) // Dynamic nodes data structure.
  {
    const Int my_idx = nodes.n(); // My node index.
    // Decide what to do.
    if (es.n() == 0) {
      // I have no elements, so return 0 to indicate I'm a leaf node containing
      // nothing.
      return 0;
    } else if (es.n() <= o.max_nelem || depth == max_depth_) {
      // I'm a leaf node with elements. Store my list of elements and return the
      // storage location.
      const Int os = offsets.back();
      offsets.push(os + es.n());
      for (Int i = 0, n = es.n(); i < n; ++i)
        elems[os + i] = es[i];
      return 1 - offsets.n();
    } else {
      // I'm not a leaf node.
      nodes(0, my_idx) = 0; // Insert myself Into the nodes array.
      for (Int ic = 0; ic < 8; ++ic) {
        BoundingBox child_bb;
        fill_child_bb(nbb, ic, child_bb);
        // Find the elements that are in this child's bb.
        IntList ces(wrk.data());
        for (Int i = 0, n = es.n(); i < n; ++i)
          if (do_bb_overlap(child_bb, slice(ebbs, es[i])))
            ces.push(es[i]);
        // Create some work space.
        IntList cwrk(wrk.data() + ces.n());
        // Recurse.
        const Int child_idx = init_r(depth+1, child_bb, ebbs, o, ces, cwrk,
                                     offsets, elems, nodes);
        nodes(ic, my_idx) = child_idx;
      }
      return my_idx;
    }
  }

  void init_static_ds (const DynNodes nodes, const DynIntList& offsets,
                       const DynIntList& elems) {
    {
      ko::resize(nodes_, nodes.n(), 8);
      auto nodes_hm = ko::create_mirror_view(nodes_);
      for (Int i = 0; i < nodes.n(); ++i)
        for (Int j = 0; j < 8; ++j)
          nodes_hm(i,j) = nodes(j,i);
      ko::deep_copy(nodes_, nodes_hm);
    }
    hm_resize_and_copy(offsets_, offsets, offsets.n());
    hm_resize_and_copy(elems_, elems, elems.n());
  }

  // Using parent bb p, fill child bb c, with child_idx in 0:7.
  template <typename CBB, typename BB>
  KOKKOS_INLINE_FUNCTION
  static void fill_child_bb (const CBB& p, const Int& child_idx, BB& c) {
    const Real m[] = { 0.5*(p[0] + p[3]),
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
  template <typename BB>
  KOKKOS_INLINE_FUNCTION
  static bool do_bb_overlap (const BoundingBox a, const BB b) {
    for (Int i = 0; i < 3; ++i)
      if ( ! do_lines_overlap(a[i], a[i+3], b[i], b[i+3]))
        return false;
    return true;
  }

  KOKKOS_INLINE_FUNCTION
  static bool do_lines_overlap (const Real& a1, const Real& a2,
                                const Real& b1, const Real& b2) {
    return ! (a2 < b1 || a1 > b2);
  }

  template <typename CV, typename Functor> KOKKOS_INLINE_FUNCTION
  void apply_r (const Int ni, const BoundingBox& nbb, const CV bb,
                Functor& f) const {
    for (Int i = 0; i < 8; ++i) {
      BoundingBox child_bb;
      fill_child_bb(nbb, i, child_bb);
      if ( ! do_bb_overlap(child_bb, bb)) continue;
      Int e = nodes_(ni,i);
      if (e > 0)
        apply_r(e, child_bb, bb, f);
      else if (e < 0) {
        e = std::abs(e + 1);
        for (Int k = offsets_[e]; k < offsets_[e+1]; ++k)
          f(elems_[k]);
      }
    }
  }
};

namespace test {
static constexpr Int max_nvert = 20;
static constexpr Int max_hits = 25; // Covers at least a 2-halo.

// In practice, we want to form high-quality normals using information about the
// mesh.
template <typename geo>
void fill_normals (sh::Mesh<ko::HostSpace>& m) {
  // Count number of edges.
  Int ne = 0;
  for (Int ip = 0; ip < nslices(m.e); ++ip)
    for (Int iv = 0; iv < szslice(m.e); ++iv)
      if (m.e(ip,iv) == -1) break; else ++ne;
  // Fill.
  Idxs::HostMirror en("en", nslices(m.e), szslice(m.e));
  ko::deep_copy(en, -1);
  Vec3s::HostMirror nml("nml", ne, 3);
  Int ie = 0;
  for (Int ip = 0; ip < nslices(m.e); ++ip)
    for (Int iv = 0; iv < szslice(m.e); ++iv)
      if (m.e(ip,iv) == -1)
        break;
      else {
        // Somewhat complicated next node index.
        const Int iv_next = (iv+1 == szslice(m.e) ? 0 :
                             (m.e(ip,iv+1) == -1 ? 0 : iv+1));
        geo::edge_normal(slice(m.p, m.e(ip, iv)), slice(m.p, m.e(ip, iv_next)),
                         slice(nml, ie));
        en(ip,iv) = ie;
        ++ie;
      }
  m.en = en;
  m.nml = nml;
}

//todo The current approach is to do redundant clips so that the hits buffer can
// be small and static. Need to think about this.
template <typename geo>
class AreaOTFunctor {
  const sh::Mesh<>& cm_;
  const ConstVec3s& p_;
  const ConstIdxs& e_;
  const Int k_; // Index into (p,e).
  //todo More efficient method that also works on GPU.
  Int hits_[max_hits];
  Int nh_;
  Real area_;

public:
  KOKKOS_INLINE_FUNCTION
  AreaOTFunctor (const sh::Mesh<>& cm, const ConstVec3s& p, const ConstIdxs& e,
                 const Int& k)
    : cm_(cm), p_(p), e_(e), k_(k), nh_(0), area_(0)
  {}

  KOKKOS_INLINE_FUNCTION void operator() (const Int mesh_elem_idx) {
    // Check whether we've clipped against this polygon before and there was a
    // non-0 intersection.
    for (Int i = 0; i < nh_; ++i)
      if (hits_[i] == mesh_elem_idx)
        return;
    // We have not, so do the intersection.
    Int no = 0;
    {
      // Area of all overlapping regions.
      // In and out vertex lists.
      Real buf[9*max_nvert];
      RawVec3s
        vi(buf, max_nvert, 3),
        vo(buf + 3*max_nvert, max_nvert, 3),
        wrk(buf + 6*max_nvert, max_nvert, 3);
      Int ni;
      ni = 0;
      for (Int i = 0; i < szslice(e_); ++i) {
        if (e_(k_,i) == -1) break;
        copy(slice(vi, i), slice(p_, e_(k_,i)), 3);
        ++ni;
      }
      sh::clip_against_poly<geo>(cm_, mesh_elem_idx, vi, ni, vo, no, wrk);
      if (no) area_ += geo::calc_area(vo, no);
    }
    if (no) {
      // Non-0 intersection, so record.
      if (nh_ == max_hits) Kokkos::abort("max_hits is too small.");
      hits_[nh_++] = mesh_elem_idx;
    }
  }

  KOKKOS_INLINE_FUNCTION const Real& area () const { return area_; }
};

template <typename geo, Int ot_max_depth_>
class TestAreaOTFunctor {
  typedef Octree<ot_max_depth_> OctreeT;

  const sh::Mesh<> cm_;
  const OctreeT ot_;
  mutable ConstVec3s p_;
  mutable ConstIdxs e_;

public:
  typedef Real value_type;

  TestAreaOTFunctor (const sh::Mesh<ko::HostSpace>& cm,
                     const ConstVec3s::HostMirror& p_hm,
                     const ConstIdxs::HostMirror& e_hm, const OctreeT& ot)
    : cm_(cm), ot_(ot)
  {
    { Vec3s p; resize_and_copy(p, p_hm); p_ = p; }
    { Idxs e; resize_and_copy(e, e_hm); e_ = e; }
  }
  
  // Clip the k'th polygon in (p,e) against mesh cm.
  KOKKOS_INLINE_FUNCTION void operator() (const Int k, Real& area) const {
    // Clipped element bounding box.
    Real ebb[6];
    OctreeT::calc_bb(p_, slice(e_, k), szslice(e_), ebb);
    // Get list of possible overlaps.
    AreaOTFunctor<geo> f(cm_, p_, e_, k);
    //todo Team threads.
    ot_.apply(ebb, f);
    area += f.area();
  }
};

template <typename geo> Real test_area_ot (
  const ConstVec3s::HostMirror& cp, const ConstIdxs::HostMirror& ce,
  const ConstVec3s::HostMirror& p, const ConstIdxs::HostMirror& e)
{
  typedef Octree<10> OctreeT;

  // Clip mesh and edge normal calculation. (In practice, we'd like to use
  // higher-quality edge normals.)
  sh::Mesh<ko::HostSpace> cm; cm.p = cp; cm.e = ce;
  fill_normals<geo>(cm);

  Real et[2];
  auto t = tic();
  // Oct-tree over the clip mesh.
  OctreeT ot(cp, ce);
  et[0] = toc(t);

  Real area = 0;
  TestAreaOTFunctor<geo, OctreeT::max_depth> f(cm, p, e, ot);
  t = tic();
  ko::parallel_reduce(nslices(e), f, area);
  et[1] = toc(t);
#ifdef SIQK_TIME
  printf("%10d", nslices(ce));
  print_times("test_area_ot", et, 2);
#endif
  return area;
}
} // namespace test
} // namespace siqk

#endif // INCLUDE_SIQK_HPP
