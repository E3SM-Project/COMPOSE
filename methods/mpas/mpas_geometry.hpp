#ifndef INCLUDE_MPAS_GEOMETRY_HPP
#define INCLUDE_MPAS_GEOMETRY_HPP

#include "mpas_types.hpp"
#include "mpas_util.hpp"
#include "mpas_triquad.hpp"

// Mostly from siqk_geometry.hpp. Where relevant, the sphere is
// assumed to have unit radius.

template <int n>
struct vec {
  template <typename V>
  static inline void set (V v, const Real& a) {
    for (int i = 0; i < n; ++i) v[i] = a;
  }
  template <typename V>
  static inline void scale (const Real& a, V v) {
    for (int i = 0; i < n; ++i) v[i] *= a;
  }
  template <typename CV, typename V>
  static inline void axpy (const Real& a, const CV x, V y) {
    for (int i = 0; i < n; ++i) y[i] += a*x[i];
  }
  template <typename CVX, typename CVY, typename VZ>
  static inline void axpbyz (const Real& a, const CVX x, const Real& b, const CVY y,
                             VZ z) {
    for (int i = 0; i < n; ++i) z[i] = a*x[i] + b*y[i];
  }
};

struct geo {
template <typename CVA, typename CVB, typename V>
static inline void cross (const CVA a, const CVB b, V c) {
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}
template <typename CVA, typename CVB>
static inline Real dot (const CVA a, const CVB b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
template <typename CV>
static inline Real norm2 (const CV v) {
  return dot(v, v);
}
template <typename CV>
static inline Real norm (const CV v) {
  return std::sqrt(dot(v, v));
}
template <typename CVA, typename CVB>
static inline Real equal (const CVA a, const CVB b) {
  return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}
template <typename V>
static inline void set (V v, const Real& a) {
  v[0] = a; v[1] = a; v[2] = a;
}
template <typename V>
static inline void setn (const int n, V v, const Real& a) {
  for (int i = 0; i < n; ++i) v[i] = 0;
}
template <typename V>
static inline void scale (const Real& a, V v) {
  v[0] *= a; v[1] *= a; v[2] *= a;
}
template <typename V>
static inline void scalen (const int n, const Real& a, V v) {
  for (int i = 0; i < n; ++i) v[i] *= a;
}
template <typename V>
static inline void normalize (V v) {
  scale(1.0/std::sqrt(norm2(v)), v);
}
template <typename CVC, typename CVA, typename CVB>
static inline Real dot_c_amb (const CVC c, const CVA a, const CVB b) {
  return c[0]*(a[0] - b[0]) + c[1]*(a[1] - b[1]) + c[2]*(a[2] - b[2]);
}
template <typename CV, typename V>
static inline void axpy (const Real& a, const CV x, V y) {
  vec<3>::axpy(a, x, y);
}
template <typename CVX, typename CVY, typename VZ>
static inline void axpbyz (const Real& a, const CVX x, const Real& b, const CVY y,
                    VZ z) {
  vec<3>::axpbyz(a, x, b, y, z);
}
template <typename CV, typename V>
static inline void axpyn (const int n, const Real& a, const CV x, V y) {
  if (n == 3) vec<3>::axpy(a, x, y);
  else if (n == 4) vec<4>::axpy(a, x, y);
  else for (int i = 0; i < n; ++i) y[i] += a*x[i];
}
template <typename CVX, typename CVY, typename VZ>
static inline void axpbyzn (const int n, const Real& a, const CVX x, const Real& b,
                            const CVY y, VZ z) {
  if (n == 3) vec<3>::axpbyz(a, x, b, y, z);
  else if (n == 4) vec<4>::axpbyz(a, x, b, y, z);
  else for (int i = 0; i < n; ++i) z[i] = a*x[i] + b*y[i];
}
template <typename V, typename CV>
static inline void copy (V d, const CV s) {
  d[0] = s[0];
  d[1] = s[1];
  d[2] = s[2];
}
template <typename V, typename CV>
static inline void copyn (const int n, V d, const CV s) {
  for (int i = 0; i < n; ++i) d[i] = s[i];
}
template <typename CV, typename V>
static inline void combine (const CV u, const CV v, const Real& a, V x) {
  const Real& oma = 1 - a;
  x[0] = oma*u[0] + a*v[0];
  x[1] = oma*u[1] + a*v[1];
  x[2] = oma*u[2] + a*v[2];
}

template <typename CV, typename V>
static inline void edge_normal (const CV a, const CV b, V en) {
  cross(a, b, en);
  normalize(en);
}

// Is v inside the line anchored at a with inward-facing normal n?
template <typename CVV, typename CVA, typename CVN>
static inline bool inside (const CVV v, const CVA a, const CVN n) {
  return dot_c_amb(n, v, a) >= 0;
}

template <typename CV, typename V>
static inline void intersect (const CV v1, const CV v2,
                              const CV e1, const CV en,
                              V intersection) {
  Real a; {
    const Real
      num = dot_c_amb(en, e1, v1),
      den = dot_c_amb(en, v2, v1);
    a = num == 0 || den == 0 ? 0 : num/den;
  }
  combine(v1, v2, a, intersection);
  normalize(intersection);
}

template <typename CV, typename V>
static bool output (const CV v, int& no, V vo, const int vo_cap) {
  if (no >= vo_cap) return false;
  const auto vo_no = slice(vo, no);
  vo_no[0] = v[0];
  vo_no[1] = v[1];
  vo_no[2] = v[2];
  ++no;
  return true;
}

template <typename CV, typename V>
static inline void intersect_from_2linesegs (const CV l1v1, const CV l1v2,
                                             const CV l2v1, const CV l2v2,
                                             V intersection) {
  Real en[3];
  edge_normal(l2v1, l2v2, en);
  intersect(l1v1, l1v2, l2v1, en, intersection);
}

// Eq after eq 10 in Lauritzen et al test cases paper.
template <typename CV>
static inline Real calc_arc_length (const CV a, const CV b) {
  const Real d = dot(a, b);
  if (d >= 1) return 0;
  Real cp[3];
  cross(a, b, cp);
  const Real cpn = std::sqrt(norm2(cp));
  return std::atan2(cpn, d);
}

// This uses a terrible formula, but it's just for testing.
template <typename CV3s>
static inline Real calc_area_formula (const CV3s& v, const int n) {
  Real area = 0;
  for (int i = 1, ilim = n - 1; i < ilim; ++i) {
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

template <typename CV3s>
static Real calc_area (const TriangleQuadrature& q, const CV3s& v,
                       const Int n, const Int order=8) {
  const Real* coord;
  const Real* weight;
  int nqp;
  q.get_coef(order, coord, weight, nqp);
  Real area = 0, u[3];
  for (Int i = 1, ilim = n - 1; i < ilim; ++i) {
    for (Int k = 0; k < nqp; ++k) {
      const Real jac = calc_tri_jacobian(slice(v,0), slice(v,i), slice(v,i+1),
                                         slice(coord, k), u, false);
      area += weight[k]*jac;
    }
  }
  return 0.5*area;
}

template <typename CV, typename CA>
static Real calc_tri_jacobian (const CV v1, const CV v2, const CV v3,
                               const CA alpha, Real u[3],
                               // if false, return < 0 for CW triangle, > 0 for CCW
                               const bool nosign = true) {
  // V(:,i) is vertex i of the spherical triangle on the unit sphere. The
  // coefs
  //     a = [alpha1, alpha2, 1 - alpha1 - alpha2]'
  //       = [1 0; 0 1; -1 -1] [alpha1, alpha2]'
  //       = a_alpha alpha
  // (barycentric coords) give the location
  //     v = V a = V a_alpha alpha
  // on the planar triangle, and u = uvec(v) is the point on the unit sphere.
  //   For a planar tri in 3D, the jacobian is
  //     v_alpha = v_a a_alpha
  //             = V [1 0; 0 1; -1 -1]
  //     J = norm(cross(v_alpha(:,1), v_alpha(:,2))).
  // For a spherical tri with the same vertices,
  //     u = v/(v' v)^{1/2}
  //     u_alpha = u_v v_alpha
  //             = (v'v)^{-1/2} (I - u u') v_alpha
  //     J = norm(cross(u_alpha(:,1), u_alpha(:,2))).
  for (Int k = 0; k < 3; ++k) u[k] = 0;
  axpy(alpha[0], v1, u);
  axpy(alpha[1], v2, u);
  axpy(alpha[2], v3, u);
  const auto oovn = 1/std::sqrt(norm2(u));
  scale(oovn, u);
  Real u_a[3][3];
  axpbyz(1, v1, -1, v3, u_a[0]);
  axpbyz(1, v2, -1, v3, u_a[1]);
  for (int i = 0; i < 2; ++i) {
    axpy(-dot(u, u_a[i]), u, u_a[i]);
    scale(oovn, u_a[i]);
  }
  cross(u_a[0], u_a[1], u_a[2]);
  const Real jac = std::sqrt(norm2(u_a[2]));
  if (nosign) return jac;
  return (dot(u_a[2], u) >= 0 ? 1 : -1)*jac;
}

template <typename CA, typename CB>
static Real reldif (const CA a, const CB b) {
  Real u[3];
  geo::axpbyz(1, b, -1, a, u);
  return std::sqrt(geo::norm2(u)/geo::norm2(a));
}
}; // struct geo

// Input is in radians.
static inline void lonlat2xyz (const Real lon, const Real lat,
                               Real& x, Real& y, Real& z,
                               const Real radius = 1) {
  const Real sinl = std::sin(lat), cosl = std::cos(lat);
  x = radius*std::cos(lon)*cosl;
  y = radius*std::sin(lon)*cosl;
  z = radius*sinl;
}

static inline void xyz2lonlat (const Real x, const Real y, const Real z,
                               Real& lon, Real& lat) {
  lon = std::atan2(y, x);
  const Real r = std::sqrt(square(x) + square(y) + square(z));
  lat = std::asin(z/r);
}

// Eq after eq 10 in Lauritzen et al test cases paper.
static inline Real great_circle_dist (
  const Real lon1, const Real lat1, const Real lon2, const Real lat2,
  const Real R = 1)
{
  Real xA, yA, zA;
  Real xB, yB, zB;
  lonlat2xyz(lon1, lat1, xA, yA, zA);
  lonlat2xyz(lon2, lat2, xB, yB, zB);
  Real cp1, cp2, cp3, cpnorm, dotprod;
  cp1 = yA*zB - yB*zA;
  cp2 = xB*zA - xA*zB;
  cp3 = xA*yB - xB*yA;
  cpnorm = std::sqrt(cp1*cp1 + cp2*cp2 + cp3*cp3);
  dotprod = xA*xB + yA*yB + zA*zB;
  return R*std::atan2(cpnorm, dotprod);
}

// Convert from (u,v), where u is velocity along latitude and v is velocity
// along longitude, to (x,y,z), which is velocity in the global cartesian
// coordinate system. Add a w (local vertical) component to push the position
// (X,Y,Z) back to the unit sphere.
inline void uv2xyz (
  const Real X, const Real Y, const Real Z, // position
  const Real u, const Real v, // velocity in tangent coord system
  Real& x, Real& y, Real& z)  // velocity in global coord system
{
  // r should be 1 but will numerically drift, so measure it ...
  const Real r = std::sqrt(X*X + Y*Y + Z*Z);
  // ... and then add a local vertical velocity to project back to the sphere.
  const Real w = 1 - r;
  Real R[9]; // Row major.
  // The local vertical is just the position vector.
  R[2] = X/r; R[5] = Y/r; R[8] = Z/r;
  // The local along-latitude vector.
  R[0] = -Y; R[3] = X; R[6] = 0;
  const Real den = std::sqrt(R[0]*R[0] + R[3]*R[3]);
  if (den != 0) { R[0] /= den; R[3] /= den; }
  // Local vertical x along-latitude.
  R[1] = R[5]*R[6] - R[8]*R[3];
  R[4] = R[8]*R[0] - R[2]*R[6];
  R[7] = R[2]*R[3] - R[5]*R[0];
  // Transform.
  x = R[0]*u + R[1]*v + R[2]*w;
  y = R[3]*u + R[4]*v + R[5]*w;
  z = R[6]*u + R[7]*v + R[8]*w;
}

#endif
