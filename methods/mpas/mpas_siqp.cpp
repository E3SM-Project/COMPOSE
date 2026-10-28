// Derived from SIQK package of https://github.com/E3SM-Project/COMPOSE
// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#include "mpas_siqp.hpp"

#include <set>
#include <iostream>
#include <limits>

namespace siqp {

#define pr(m) do {                              \
    std::stringstream _ss_;                     \
    _ss_ << m << std::endl;                     \
    std::cerr << _ss_.str();                    \
  } while (0)
#define prc(m) pr(#m << " | " << (m))
#define puf(m)"(" << #m << " " << (m) << ")"
#define pu(m) << " " << puf(m)
template<typename T>
static void prarr (const std::string& name, const T* const v, const size_t n) {
  std::cerr << name << ": ";
  for (size_t i = 0; i < n; ++i) std::cerr << " " << v[i];
  std::cerr << "\n";
}

#ifdef SIQP_TIME
timeval tic () {
  timeval t;
  gettimeofday(&t, 0);
  return t;
}
Real calc_et (const timeval& t1, const timeval& t2) {
  static const Real us = 1.0e6;
  return (t2.tv_sec * us + t2.tv_usec - t1.tv_sec * us - t1.tv_usec) / us;
}
Real toc (const timeval& t1) {
  timeval t;
  gettimeofday(&t, 0);
  return calc_et(t1, t);
}
Real get_memusage () {
  static const Real scale = 1.0 / (1 << 10); // Memory in MB.
  rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  return ru.ru_maxrss*scale;
}
#endif
void print_times (const std::string& name, const Real* const parts,
                         const int nparts) {
#ifdef SIQP_TIME
  Real total = 0; for (int i = 0; i < nparts; ++i) total += parts[i];
  printf("%20s %1.3e s %7.1f MB", name.c_str(), total, get_memusage());
  for (int i = 0; i < nparts; ++i) printf(" %1.3e s", parts[i]);
  printf("\n");
#endif
}

int solve_quadratic_equation (const Real& a, const Real& b, const Real& c,
                              Real xs[2]) {
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

// Geometry on the sphere. All inputs and outputs are relative to the
// unit-radius sphere.
struct SphereGeometry {
  enum { dim = 3 };

  template <typename CV, typename V>
  static void cross (const CV a, const CV b, V c) {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
  }
  template <typename CA, typename CB>
  static Real dot (const CA a, const CB b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  }
  template <typename CV>
  static Real norm2 (const CV v) {
    return dot(v, v);
  }
  template <typename V>
  static void scale (const Real a, V v) {
    v[0] *= a; v[1] *= a; v[2] *= a;
  }
  template <typename V>
  static void set (V v, const Real a) {
    v[0] = a; v[1] = a; v[2] = a;
  }
  template <typename V>
  static void normalize (V v) {
    scale(1.0/std::sqrt(norm2(v)), v);
  }
  template <typename CV>
  static Real dot_c_amb (const CV c, const CV a, const CV b) {
    return c[0]*(a[0] - b[0]) + c[1]*(a[1] - b[1]) + c[2]*(a[2] - b[2]);
  }
  template <typename CV, typename V>
  static inline void axpy (const Real& a, const CV x, V y) {
    y[0] += a*x[0];
    y[1] += a*x[1];
    y[2] += a*x[2];
  }
  template <typename CVX, typename CVY, typename VZ>
  static inline void axpbyz (const Real& a, const CVX x, const Real& b, const CVY y,
                             VZ z) {
    z[0] = a*x[0] + b*y[0];
    z[1] = a*x[1] + b*y[1];
    z[2] = a*x[2] + b*y[2];
  }
  template <typename V, typename CV>
  static void copy (V d, const CV s) {
    d[0] = s[0];
    d[1] = s[1];
    d[2] = s[2];
  }
  template <typename CV, typename V>
  static void combine (const CV u, const CV v, const Real a, V x) {
    const Real oma = 1 - a;
    x[0] = oma*u[0] + a*v[0];
    x[1] = oma*u[1] + a*v[1];
    x[2] = oma*u[2] + a*v[2];
  }

  template <typename CV, typename V>
  static void edge_normal (const CV a, const CV b, V en) {
    cross(a, b, en);
    normalize(en);
  }

  // Is v inside the line (a1,a2) having normal n?
  template <typename CV>
  static bool inside (const CV v, const CV a1, const CV n) {
    return dot_c_amb(n, v, a1) > 0;
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
  template <typename CV, typename V>
  static void intersect (const CV v1, const CV v2, const CV e1, const CV en,
                         V intersection) {
    Real a; {
      const Real
        num = dot_c_amb(en, e1, v1),
        den = dot_c_amb(en, v2, v1);
      a = num == 0 || den == 0 ? 0 : num/den;
      a = a < 0 ? 0 : a > 1 ? 1 : a;
    }
    combine(v1, v2, a, intersection);
    normalize(intersection);
  }

  template <typename CV>
  static bool output (const CV v, int& no, Array2D<Real>& vo) {
#ifndef NDEBUG
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

  // This uses a terrible formula, but it's just for testing. 
  static Real calc_area (const Array2D<const Real>& v) {
    Real area = 0;
    for (int i = 1; i < v.n() - 1; ++i) {
      const Real a = calc_arc_length(v(0), v(i));
      const Real b = calc_arc_length(v(i), v(i+1));
      const Real c = calc_arc_length(v(i+1), v(0));
      const Real s = 0.5*(a + b + c);
      const Real d = (std::tan(0.5*s)*std::tan(0.5*(s-a))*
                      std::tan(0.5*(s-b))*std::tan(0.5*(s-c)));
      if (d <= 0) continue;
      area += 4*std::atan(std::sqrt(d));
    }
    return area;
  }
  template <typename CV>
  static Real calc_arc_length (const CV a, const CV b) {
    const Real d = dot(a, b);
    if (d >= 1) return 0;
    return std::acos(d);
  }

  /* For quadratic edges. Let
       n = edge normal
       c = edge point
       d = n'c
       v(a) = (1-a)^2 s + a (1-a) m + a^2 p
     Solve n' v = d for a:
       (1-a)^2 s'n + a (1-a) m'n + a^2 p'n = n'c
       (s - m + p)'n a^2 + (-2 s + m)'n a + (s - c)'n = 0
  */
  template <typename CV, typename V>
  static int intersect (const CV s, const CV m, const CV p, const CV e1,
                        const CV nml, V as) {
    Real w[3];
    for (int i = 0; i < 3; ++i) w[i] = s[i] - m[i] + p[i];
    const Real a = dot(nml, w);
    for (int i = 0; i < 3; ++i) w[i] = m[i] - 2*s[i];
    const Real b = dot(nml, w);
    for (int i = 0; i < 3; ++i) w[i] = s[i] - e1[i];
    const Real c = dot(nml, w);
    return solve_quadratic_equation(a, b, c, as);
  }
};

typedef SphereGeometry geo;

namespace sh {
template <typename CV>
bool clip_against_edge (
  // Input vertex list.
  const Array2D<const Real>& vi, const int ni,
  // Output vertex list.
  Array2D<Real>& vo, int& no,
  // The end points of the clip edge segment.
  const CV ce1,
  // Clip edge's inward-facing normal.
  const CV cen)
{
  const Real* s, * p;
  Real intersection[3];
  no = 0;
  s = vi(ni-1);
  for (int j = 0; j < ni; ++j) {
    p = vi(j);
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

bool clip_against_poly (
  const Mesh& m, const int cp_e,
  const Array2D<const Real>& vi, const int ni,
  Array2D<Real>& vo, int& no,
  Real* const wrk, const int nvertwrk)
{
  Array2D<Real> vo1(3, nvertwrk, wrk);
  int nos[] = { 0, 0 };
  Array2D<Real>* vs[] = { &vo, &vo1 };

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

  if ( ! clip_against_edge(vi, ni, *vs[0], nos[0], m.p(e[0]), m.nml(en[0])))
    return false;
  if ( ! nos[0]) return true;

  for (int ie = 1, ielim = nv - 1; ; ++ie) {
    if ( ! clip_against_edge(*vs[0], nos[0], *vs[1], nos[1], m.p(e[ie]),
                             m.nml(en[ie])))
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
template <typename geo>
bool clip_against_poly (
  // Clip polygon.
  const Array2D<const Real>& clip_poly,
  // Clip polygon edges' inward-facing normals.
  const Array2D<const Real>& clip_edge_normals,
  const Array2D<const Real>& vi, const int ni,
  Array2D<Real>& vo, int& no,
  Real* const wrk, const int nvertwrk)
{
  Array2D<Real> vo1(3, nvertwrk, wrk);
  int nos[] = { 0, 0 };
  Array2D<Real>* vs[] = { &vo, &vo1 };

  no = 0;
  const auto nv = clip_poly.n();
  if (nv % 2 == 0) {
    // Make sure the final vertex output list is in the caller's buffer.
    std::swap(vs[0], vs[1]);
    std::swap(nos[0], nos[1]);
  }

  if ( ! clip_against_edge(vi, ni, *vs[0], nos[0], clip_poly(0),
                           clip_edge_normals(0)))
    return false;
  if ( ! nos[0]) return true;

  for (int ie = 1, ielim = nv - 1; ; ++ie) {
    if ( ! clip_against_edge(*vs[0], nos[0], *vs[1], nos[1], clip_poly(ie),
                             clip_edge_normals(ie)))
      return false;
    if ( ! nos[1]) return true;
    if (ie == ielim) break;
    std::swap(vs[0], vs[1]);
    std::swap(nos[0], nos[1]);
  }

  no = nos[1];
  return true;
}

template <typename geo>
struct ice {
  typedef double Real;
  typedef int Int;
  typedef unsigned int UInt;
  enum { dim = geo::dim };

  template <typename CA, typename CB>
  static Real reldif (const CA a, const CB b) {
    Real u[dim];
    geo::axpbyz(1, b, -1, a, u);
    return std::sqrt(geo::norm2(u)/geo::norm2(a));
  }

  // Follow the outer loop of the Sutherland-Hodgmann algorithm.
  static bool clip_cpoly_against_convex_poly (
    const Mesh& m, const int cp_e,
    const Array2D<const Real>& vi, const Array1D<const Int>& vti, const Int ni,
    Array2D<Real>& vo, Array1D<Int>& vto, Int& no,
    Real* const rwrk, Int* const iwrk, const Int nwrk)
  {
    Array2D<Real> vo1(dim, nwrk, rwrk);
    Array2D<Real>* vs[] = { &vo, &vo1 };
    Array1D<Int> vto1(nwrk, iwrk);
    Array1D<Int>* vts[] = { &vto, &vto1 };
    int nos[] = { 0, 0 };

    const auto e = m.e(cp_e);
    const auto en = m.en(cp_e);
  
    auto nv = m.e.m(); // Number of vertices in clip polygon.
    while (e[nv-1] == -1) --nv;

    no = 0;
    if (nv % 2 == 0) {
      // Make sure the final vertex output list is in the caller's buffer.
      std::swap(vs[0], vs[1]);
      std::swap(vts[0], vts[1]);
      std::swap(nos[0], nos[1]);
    }

    if ( ! clip_cpoly_against_sed(vi, vti, ni, *vs[0], *vts[0], nos[0],
                                  m.p(e[0]), m.nml(en[0])))
      return false;
    if ( ! nos[0]) return true;

    for (Int ie = 1, ielim = nv - 1; ; ++ie) {
      if ( ! clip_cpoly_against_sed(*vs[0], *vts[0], nos[0], *vs[1], *vts[1], nos[1],
                                    m.p(e[ie]), m.nml(en[ie])))
        return false;
      if ( ! nos[1]) return true;
      if (ie == ielim) break;
      std::swap(vs[0], vs[1]);
      std::swap(vts[0], vts[1]);
      std::swap(nos[0], nos[1]);
    }

    no = nos[1];
    return true;    
  }

  template <typename CV>
  static bool clip_cpoly_against_sed (
    // Input vertex and vertex type lists.
    const Array2D<const Real>& vi, const Array1D<const Int>& vti, const Int ni,
    // Outputs.
    Array2D<Real>& vo, Array1D<Int>& vto, Int& no,
    // An end point of the clip edge segment.
    const CV se1,
    // Clip edge's inward-facing normal.
    const CV sen)
  {
    bool ends_connected = false;
    Int k = 0;
    no = 0;
    for (;;) {
      if (no + 5 > nslices(vo)) return false;
      if (vti[k] == 1) {
        // 1. Start of vi list, with [1 0] as the two first vts.
        assert(k == 0);
        assert(vti[ni-1] == 0);
        no += clip_ced_against_sed(se1, sen,
                                   slice(vi, ni-1), slice(vi, 0), slice(vi, 1),
                                   vo, vto);
        k++;
        ends_connected = true;
      } else if (vti[k] == 0 && k+1 == ni) {
        // 2. Last vertex in the list.
        if (ends_connected) break;
        if (vti[0] == 0) {
          // 2a. A sedge connects end of list to start.
          const UInt nv = clip_sed_against_sed(se1, sen,
                                               slice(vi, ni-1), slice(vi, 0),
                                               offset(vo, no));
          for (UInt i = 0; i < nv; ++i) vto[no+i] = 0;
          no += nv;
        } else {
          // 2b. Block 1 took care of this cedge.
          assert(false);
        }
        break;
      } else if (vti[k+1] == 1) {
        assert(vti[k] == 0);
        if (k+2 == ni) {
          // 3. [0 1] at the end of the list connects to 0 at the start.
          assert(vti[0] == 0);
          no += clip_ced_against_sed(se1, sen,
                                     slice(vi, k), slice(vi, k+1), slice(vi, 0),
                                     offset(vo, no), offset(vto, no));
          break;
        } else {
          assert(k+2 < ni);
          // 4. General case: [0 1 0] in the middle of the list.
          no += clip_ced_against_sed(se1, sen,
                                     slice(vi, k), slice(vi, k+1), slice(vi, k+2),
                                     offset(vo, no), offset(vto, no));
          k += 2;
        }
      } else {
        // 5. General case: [0 0] in the middle of the list.
        assert(k+1 < ni);
        assert(vti[k+1] == 0);
        const UInt nv = clip_sed_against_sed(se1, sen,
                                             slice(vi, k), slice(vi, k+1),
                                             offset(vo, no));
        for (UInt i = 0; i < nv; ++i) vto[no+i] = 0;
        no += nv;
        k++;
      }
    }
    return true;
  }

  // sed must have >= 2 slices allocated.
  template <typename CV, typename Array>
  static UInt clip_sed_against_sed (const CV se1, const CV sen,
                                    const CV s, const CV p, Array sed) {
    const bool s_inside = geo::inside(s, se1, sen);
    const bool p_inside = geo::inside(p, se1, sen);
    if (p_inside) {
      if (s_inside) {
        copy(slice(sed, 0), p);
        return 1;
      } else {
        geo::intersect(s, p, se1, sen, slice(sed, 0));
        copy(slice(sed, 1), p);
        return 2;
      }
    } else {
      if (s_inside) {
        geo::intersect(s, p, se1, sen, slice(sed, 0));
        return 1;
      } else {
        return 0;
      }
    }
  }

  // eds and evts must have >= 5 slices allocated.
  template <typename CV, typename Array, typename IV>
  static UInt clip_ced_against_sed (const CV se1, const CV sen,
                                    const CV s, const CV m, const CV p,
                                    Array eds, IV edvts) {
    Real as[2] = {0};
    UInt nas = intersect(s, m, p, se1, sen, as);
    bool s_inside = geo::inside(s, se1, sen);
    bool p_inside = geo::inside(p, se1, sen);

    // Handle cases where FP fails to lead to a consistent state.
    if (p_inside != s_inside) {
      if (nas == 0) {
        // There is no FP intersection, so p and s might as well both be inside.
        p_inside = s_inside = true;
        // Now decide what to do with m.
        Real M[3];
        eval(s, m, p, 0.5, M);
        if ( ! geo::inside((const Real*) M, se1, sen)) {
          // s and p are barely inside and M is outside, so this is, up to FP, a
          // straight edge.
          copy(slice(eds, 0), p); edvts[0] = 0;
          return 1;
        } // Otherwise, continue to the first case below.
      } else if (nas == 2) {
        // Take the one closest to the middle.
        if (std::abs(as[0] - 0.5) > std::abs(as[1] - 0.5)) as[0] = as[1];
        nas = 1;
      }
    }

    if (p_inside) {
      if (s_inside) {
        if (nas < 2) {
          copy(slice(eds, 0), m); edvts[0] = 1;
          copy(slice(eds, 1), p); edvts[1] = 0;
          return 2;
        } else {
          middle_for_segment_0c(s, m, as[0], slice(eds, 0));
          eval(s, m, p, as[0], slice(eds, 1));
          eval(s, m, p, as[1], slice(eds, 2));
          middle_for_segment_d1(m, p, as[1], slice(eds, 3));
          copy(slice(eds, 4), p);
          edvts[0] = 1; edvts[1] = 0; edvts[2] = 0; edvts[3] = 1; edvts[4] = 0;
          return 5;
        }
      } else {
        assert(nas == 1);
        eval(s, m, p, as[0], slice(eds, 0));
        middle_for_segment_d1(m, p, as[0], slice(eds, 1));
        copy(slice(eds, 2), p);
        edvts[0] = 0; edvts[1] = 1; edvts[2] = 0;
        return 3;
      }
    } else {
      if (s_inside) {
        assert(nas == 1);
        middle_for_segment_0c(s, m, as[0], slice(eds, 0));
        eval(s, m, p, as[0], slice(eds, 1));
        edvts[0] = 1; edvts[1] = 0;
        return 2;
      } else {
        if (nas < 2) {
          return 0;
        } else {
          eval(s, m, p, as[0], slice(eds, 0));
          middle_for_segment(s, m, p, as[0], as[1], slice(eds, 1));
          eval(s, m, p, as[1], slice(eds, 2));
          edvts[0] = 0; edvts[1] = 1; edvts[2] = 0;
          return 3;
        }
      }
    }

    assert(0);
    return 0;
  }

  template <typename CV>
  static UInt intersect (const CV s, const CV m, const CV p,
                         const CV e1, const CV en, Real as[2]) {
    Int nas = geo::intersect(s, m, p, e1, en, as);
    if (nas == 2 && (as[1] < 0 || as[1] > 1)) --nas;
    if (nas >= 1 && (as[0] < 0 || as[0] > 1)) { as[0] = as[1]; --nas; }
    if (nas == 2 && as[0] > as[1]) std::swap(as[0], as[1]);
    assert(nas >= 0 && nas <= 2);
    return static_cast<UInt>(nas);
  }

  // Create m in (s,m,p) for the segment a in [c,d].
  template <typename CV, typename V>
  static void middle_for_segment (const CV s, const CV m, const CV p,
                                  const Real c, const Real d, V n) {
    if (c == 0) middle_for_segment_0c(s, m, d, n);
    else if (d == 1) middle_for_segment_d1(m, p, c, n);
    else middle_for_segment_cd(s, m, p, c, d, n);
  }
  template <typename CV, typename V>
  static void middle_for_segment_0c (const CV s, const CV m, const Real d,
                                     V n) {
    const Real tomd = 2*(1 - d);
    for (UInt i = 0; i < dim; ++i)
      n[i] = tomd*s[i] + d*m[i];
  }
  template <typename CV, typename V>
  static void middle_for_segment_d1 (const CV m, const CV p, const Real c,
                                     V n) {
    const Real omc = 1 - c, tc = 2*c;
    for (UInt i = 0; i < dim; ++i)
      n[i] = omc*m[i] + tc*p[i];
  }
  template <typename CV, typename V>
  static void middle_for_segment_cd (const CV s, const CV m, const CV p,
                                     const Real c, const Real d, V n) {
    const Real cd = c*d, tcd = 2*cd, c0 = 2*(cd - c - d + 1), c1 = c + d - tcd;
    for (UInt i = 0; i < dim; ++i)
      n[i] = c0*s[i] + c1*m[i] + tcd*p[i];
  }

  // Create m in (s,m,p) so that the curve hits n.
  template <typename CV, typename V>
  static void middle_matches (const CV s, const CV p, const CV n, V m) {
    for (UInt i = 0; i < dim; ++i)
      m[i] = 4*n[i] - s[i] - p[i];
  }

  template <typename CV, typename V>
  static void eval (const CV s, const CV m, const CV p, const Real a, V v) {
    const Real oma = 1 - a, oma2 = oma*oma, omaa = oma*a, a2 = a*a;
    for (UInt i = 0; i < dim; ++i)
      v[i] = oma2*s[i] + omaa*m[i] + a2*p[i];
  }

  template <typename V, typename CV>
  static void copy (V d, const CV s) {
    for (UInt i = 0; i < dim; ++i) d[i] = s[i];
  }
}; // struct ice

bool clip_cpoly_against_convex_poly (
  const Mesh& m, const int cp_e,
  const Array2D<const Real>& vi, const Array1D<const int>& vti, const int ni,
  Array2D<Real>& vo, Array1D<int>& vto, int& no,
  Real* const rwrk, int* const iwrk, const int nwrk)
{
  return ice<SphereGeometry>::clip_cpoly_against_convex_poly(
    m, cp_e, vi, vti, ni, vo, vto, no, rwrk, iwrk, nwrk);
}

void sedge_eval (const Real s[3], const Real p[3], const Real& a, Real u[3]) {
  SphereGeometry::combine(s, p, a, u);
  SphereGeometry::normalize(u);
}

void cedge_eval (const Real s[3], const Real m[3], const Real p[3], const Real& a,
                 Real u[3]) {
  ice<SphereGeometry>::eval(s, m, p, a, u);
  SphereGeometry::normalize(u);
}

void middle_matches (const Real s[3], const Real p[3], const Real n[3], Real m[3]) {
  ice<SphereGeometry>::middle_matches(s, p, n, m);
}

/* Let
       vs = [s m p v]
       vts = [0 1 0 0]
   be a canonicalized triangle, a vctri, v for <) shape. Let
       a = [a1, a2, 1 - a1 - a2]
   be the barycentric coordinates for s, p, v. Let
       b = a2/(a1 + a2)
       h(b) = parameterized curve for the curved edge
       v(a) = (1 - a1 - a2) v + (a1 + a2) h(b)
       u(a) = uvec(v(a))
   Then
       1 - b = a1/(a1 + a2)
       b_a1 = -a2/(a1 + a2)^2 = -b/(a1 + a2)
       b_a2 = 1/(a1 + a2) - a2/(a1 + a2)^2
            = 1/(a1 + a2) (1 - a2/(a1 + a2))
            = 1/(a1 + a2) (a1/(a1 + a2))
            = a1/(a1 + a2)^2
            = (1 - b)/(a1 + a2)
       (a1 + a2) b_a1 = -b
       (a1 + a2) b_a2 = 1 - b
       v_a1 = -v + h(b) + (a1 + a2) h_b b_a1
            = -v + h(b) - b h_b
       v_a2 = -v + h(b) + (1-b) h_b
            = v_a1 + h_b
       u_a = u_v v_a
           = u_v v_a
           = 1/sqrt(v'v) (I - u u') v_a
       jacobian determinant = norm(cross(u_a1, u_a2), 2).
   Now, specifically, let
       h(b) = (1-b)^2 s + b (1-b) m + b^2 p
   Then
       h_b = -2 (1-b) s + (1 - 2b) m + 2 b p
       v_a1 = -v + (1-b)^2 s + b (1-b) m + b^2 p - b (-2 (1-b) s + (1 - 2b) m + 2 b p)
            = -v + (1-b)^2 s + b (1-b) m + b^2 p + 2 b (1-b) s - b (1 - 2b) m - 2 b^2 p
            = -v + (1-b)^2 s + b (1-b) m + b^2 p + 2 b (1-b) s - b (1 - 2b) m - 2 b^2 p
            = -v + [(1-b)^2 + 2 b (1-b)] s + [b (1-b) - b (1 - 2b)] m + [b^2 - 2 b^2] p
            = -v + (1 - b) (1 + b) s + b^2 m - b^2 p
            = s - v - b^2 (s - m + p)
       v_a2 = v_a1 - 2 (1-b) s + (1 - 2b) m + 2 b p
            = v_a1 + 2 b (s - m + p) + (m - 2 s)
            = m - s - v + b (2 - b) (s - m + p)
   As a simple check, if M = (s + p)/2 and so m = 4 M - s - p = s + p, then
       v_a1 = s - v
       v_a2 = p - v,
   as in the straight-edged triangle case.
*/
Real calc_vctri_jacobian (const Real vtx[3], const Real s[3], const Real m[3], const Real p[3],
                          const Real alpha[2], Real u[3], const bool nosign) {
  typedef SphereGeometry geo;
  typedef ice<geo> ices;
  // curved edge parameterized curve parameter
  const Real
    a1_p_a2 = alpha[0] + alpha[1],
    b = a1_p_a2 <= 0 ? 0.5 : alpha[1]/a1_p_a2;
  // s - m + p
  Real smp[3];
  geo::axpbyz(1, s, -1, m, smp); geo::axpy(1, p, smp);
  // v_a1, v_a2
  Real v_a[2][3];
  geo::axpbyz(1, s, -1, vtx, v_a[0]); geo::axpy(-b*b, smp, v_a[0]);
  geo::axpbyz(1, m, -1, s, v_a[1]); geo::axpy(-1, vtx, v_a[1]); geo::axpy(b*(2-b), smp, v_a[1]);
  // u
  Real h[3];
  ices::eval(s, m, p, b, h);
  geo::axpbyz(1 - a1_p_a2, vtx, a1_p_a2, h, u);
  const Real v_nrm = std::sqrt(geo::norm2(u)), r_v_nrm = 1/v_nrm;
  geo::scale(r_v_nrm, u);
  // sphere part
  for (int i = 0; i < 2; ++i) {
    geo::axpy(-geo::dot(u, v_a[i]), u, v_a[i]);
    geo::scale(r_v_nrm, v_a[i]);
  }
  // jacobian determinant
  Real c[3];
  geo::cross(v_a[0], v_a[1], c);
  const Real jac = std::sqrt(geo::norm2(c));
  // optionally signed to handle self-intersection and anti-oriented polys
  if (nosign) return jac;
  return (geo::dot(c, u) >= 0 ? 1 : -1)*jac;
}

double urand () { return rand() / ((double) RAND_MAX + 1.0); }

bool test_ice_math () {
  typedef Real Vec[3];
  typedef SphereGeometry geo;
  typedef ice<SphereGeometry> ices;

  Vec s, M, m, p;
  for (int i = 0; i < 3; ++i) s[i] = urand();
  for (int i = 0; i < 3; ++i) M[i] = urand();
  for (int i = 0; i < 3; ++i) p[i] = urand();
  geo::normalize(s); geo::normalize(M); geo::normalize(p);
  ices::middle_matches(s, p, M, m);

  const Real tol = 10*std::numeric_limits<Real>::epsilon();
  const Real c = 0.42, d = 0.7;
  Vec r, q, n;
  ices::eval(s, m, p, c, r);
  ices::eval(s, m, p, d, q);
  ices::middle_for_segment(s, m, p, c, d, n);
  bool ok = true;
  for (const Real b : {0.0, 0.2, 0.3, 0.7, 1.0}) {
    const Real a = c + (d - c)*b;
    Vec va, vb;
    ices::eval(s, m, p, a, va); geo::normalize(va);
    ices::eval(r, n, q, b, vb); geo::normalize(vb);
    if (ices::reldif(va, vb) >= tol) ok = false;
  }

  Vec n1, n2;
  ices::middle_for_segment_0c(s, m, 0.3, n1);
  ices::middle_for_segment(s, m, p, 0, 0.3, n2);
  ok = ok && ices::reldif(n1, n2) < tol;
  ices::middle_for_segment_d1(m, p, 0.7, n1);
  ices::middle_for_segment(s, m, p, 0.7, 1, n2);
  ok = ok && ices::reldif(n1, n2) < tol;
  
  return ok;
}
} // namespace sh

typedef Octree::BoundingBox BoundingBox;

Octree::Octree (const Array2D<const Real>& p, const Array2D<const int>& e,
                const Options& o) {
  init(p, e, o);
}

Octree::Octree (const Array2D<const Real>& p, const Array2D<const int>& e) {
  Options o;
  init(p, e, o);
}

void Octree::apply (const BoundingBox& bb, Functor& f) const {
  if (nodes_.n() == 0) {
    for (int i = 0; i < offset_[1]; ++i)
      f(elems_[i]);
    return;
  }
  apply_r(0, bb_, bb, f);
}

// Using parent bb p, fill child bb c, with child_idx in 0:7.
static void fill_child_bb (const BoundingBox& p, const int& child_idx,
                           BoundingBox& c) {
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

static bool do_lines_overlap (const Real& a1, const Real& a2,
                              const Real& b1, const Real& b2) {
  return ! (a2 < b1 || a1 > b2);
}

// Do bounding boxes a and b overlap?
static bool do_bb_overlap (const BoundingBox a, const BoundingBox b) {
  for (int i = 0; i < 3; ++i)
    if ( ! do_lines_overlap(a[i], a[i+3], b[i], b[i+3]))
      return false;
  return true;
}

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

static int init_r (const int depth, // Tree's depth at this point, including root.
                   const BoundingBox& nbb, // My bounding box.
                   const Array2D<const Real>& ebbs, // All elements' bounding boxes.
                   const Octree::Options& o, // Options controlling construct of the tree.
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

void Octree::init (const Array2D<const Real>& p, const Array2D<const int>& e,
                   const Options& o) {
  if (e.n() == 0) return;
  // Get OT's bounding box.
  calc_bb(p, bb_);
  // Get elements' bounding boxes.
  Array2D<Real> ebbs(6, e.n());
  calc_bb(p, e, ebbs);
  // Static element lists for work. Each level has active work space.
  std::vector<int> buf(o.max_depth*e.n());
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
  // ...
}

void Octree::apply_r (const int ni, const BoundingBox& nbb, const BoundingBox& bb,
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

// In practice, we want to form high-quality normals using information about the
// mesh, such as that it is a CS mesh. For testing, form the normals from edge
// vertices. (This leads to increasing cancellation error with mesh refinement.)
void fill_normals (sh::Mesh& m) {
  // Count number of edges.
  int ne = 0;
  for (auto ip = zero(m.e); ip < m.e.n(); ++ip)
    for (auto iv = zero(m.e); iv < m.e.m(); ++iv)
      if (m.e(iv,ip) == -1) break; else ++ne;
  // Fill.
  Array2D<int> en(m.e.m(), m.e.n());
  en.set(-1);
  Array2D<Real> nml(3, ne);
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

namespace test {
static constexpr int max_nvert = 20;

// Used in Octree::apply to gather a set of possibly intersecting polygons.
struct OTSearchFunctor : public Octree::Functor {
  std::set<int> hits;
  void operator() (const int i) override { hits.insert(i); }
};

// Find the area of the overlapping part of two meshes by summing over the areas
// of the common refinement polygons. Obviously a silly thing to do, but a good
// test and demonstration problem.
class TestAreaOTFunctor {
  // Mesh against which to clip. ("Eulerian mesh".)
  sh::Mesh cm;
  // Mesh of clipped polygons. ("Departure mesh".)
  const Array2D<const Real> p; // 3 x #verts array of polygon vertices.
  const Array2D<const int> e;    // Array of polygons. e(:,k) is the k'th polygon.
  // Already initialized octree used to search for possibly intersecting
  // polygons.
  Octree ot;
  AreaFunctor* af;

public:
  typedef Real value_type;

  TestAreaOTFunctor (const sh::Mesh& cm_, const Array2D<const Real>& p_,
                     const Array2D<const int>& e_, const Octree& ot_,
                     AreaFunctor* af_)
    : cm(cm_), p(p_), e(e_), ot(ot_), af(af_)
  {}
  
  // k indexes (p,e).
  void operator() (const int k, Real& area) const {
    // Clipped element bounding box.
    Real ebb[6];
    Octree::calc_bb(p, e(k), e.m(), ebb);
    // Get list of possible overlaps.
    OTSearchFunctor f;
    ot.apply(ebb, f);
    // In and out vertex lists.
    Real buf[6*max_nvert];
    Array2D<Real>
      vi(3, max_nvert, buf),
      vo(3, max_nvert, buf + 3*max_nvert);
    int ni, no;
    // Create the polygon to be clipped.
    ni = 0;
    for (int i = 0; i < e.m(); ++i) {
      if (e(i,k) == -1) break;
      copy(vi(i), p(e(i,k)), 3);
      ++ni;
    }
    // Workspace.
    Real wrk[3*max_nvert];
    // Area of all overlapping regions.
    Real a = 0;
    for (const auto icp : f.hits) {
      sh::clip_against_poly(cm, icp, vi, ni, vo, no, wrk, max_nvert);
      if (no) {
        // A non-0 intersection was found. Accumulate the area.
        const auto this_area = (af ?
                                (*af)(vo.data(), no) :
                                geo::calc_area(Array2D<const Real>(vo.m(), no, vo.data())));
        a += this_area;
      }
    }
    // Add our area to the reduction.
    area += a;
  }
};

class IceTestAreaOTFunctor {
  sh::Mesh cm;
  const Array2D<const Real> p;
  const Array2D<const int> e;
  Octree ot;
  AreaFunctor* af;
  bool verbose;

public:
  typedef Real value_type;

  IceTestAreaOTFunctor (const sh::Mesh& cm_, const Array2D<const Real>& p_,
                        const Array2D<const int>& e_, const Octree& ot_,
                        AreaFunctor* af_, const bool verbose_)
    : cm(cm_), p(p_), e(e_), ot(ot_), af(af_), verbose(verbose_)
  {}
  
  // k indexes (p,e).
  void operator() (const int k, Real& area) const {
    // Clipped element bounding box.
    Real ebb[6];
    Octree::calc_bb(p, e(k), e.m(), ebb);
    // Get list of possible overlaps.
    OTSearchFunctor f;
    ot.apply(ebb, f);
    int ni, no;
    // Area of all overlapping regions.
    Real a = 0;
    for (const auto icp : f.hits) {
      ni = 0;
      static const int N = 2*max_nvert;
      Real rbuf[9*N];
      Array2D<Real> vi(3, N, rbuf), ivo(3, N, rbuf + 3*N), vo(3, N, rbuf + 6*N);
      int ibuf[2*N];
      Array1D<int> vti(N, ibuf), vto(N, ibuf + N);
      for (int i = 0; i < e.m(); ++i) {
        if (e(i,k) == -1) break;
        vti[2*i] = 0; vti[2*i+1] = 1;
        geo::copy(vi(2*i), p(e(i,k)));
        ni += 2;
      }
      for (int i = 0; i < e.m(); ++i) {
        Real n[3];
        geo::combine(vi(2*i), vi((2*(i+1)) % ni), 0.5, n);
        sh::ice<SphereGeometry>::middle_matches(vi(2*i), vi((2*(i+1)) % ni), n, vi(2*i+1));
      }
      Real rwrk[3*N];
      int iwrk[N];
      sh::clip_cpoly_against_convex_poly(cm, icp, vi, vti, ni, ivo, vto, no,
                                         rwrk, iwrk, N);
      int n = 0;
      for (int i = 0; i < no; ++i)
        if (true || vto[i] == 0) {
          geo::copy(vo(n), ivo(i));
          geo::normalize(vo(n));
          ++n;
        }
      no = n;
      if (no) {
        //todo Need to make a curved-edge version of this using vto.
        const Real a1 = (af ?
                         (*af)(vo.data(), no) :
                         geo::calc_area(Array2D<const Real>(vo.m(), no, vo.data())));
        a += a1;
      }
    }
    // Add our area to the reduction.
    area += a;
  }
};

Real test_area_ot (
  const Array2D<const Real>& cp, const Array2D<const int>& ce,
  const Array2D<const Real>& p, const Array2D<const int>& e,
  const bool verbose, const bool use_ice, AreaFunctor* af)
{
  // Clip mesh and edge normal calculation. (In practice, we'd like to use
  // higher-quality edge normals.)
  sh::Mesh cm; cm.p = cp; cm.e = ce;
  fill_normals(cm);

  Real et[2];
  auto t = tic();
  // Build an octree over the clip mesh.
  Octree ot(cp, ce);
  et[0] = toc(t);

  // Compute the area in a silly way to test search and interesection.
  t = tic();
  Real area = 0;
  if (use_ice)
    parallel_reduce(e.n(), test::IceTestAreaOTFunctor(cm, p, e, ot, af, verbose), area);
  else
    parallel_reduce(e.n(), test::TestAreaOTFunctor(cm, p, e, ot, af), area);
  et[1] = toc(t);
  if (verbose)
    print_times("test_area_ot", et, 2);
  return area;
}

} // namespace test
} // namespace siqp
