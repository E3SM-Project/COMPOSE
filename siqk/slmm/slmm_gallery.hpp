#ifndef INCLUDE_SLMM_GALLERY_HPP
#define INCLUDE_SLMM_GALLERY_HPP

#include "slmm_defs.hpp"
#include "slmm_time_int.hpp"

namespace slmm {
namespace gallery {

class OdeFnBasicRecorder {
  mutable int ne_;
  bool xyz_form_;
public:
  OdeFnBasicRecorder () : ne_(0), xyz_form_(false) {}
  void record (const Real t, const Real* const y) const { ++ne_; }
  int ne () const { return ne_; }
  void set_xyz_form (const bool use_xyz_form) { xyz_form_ = use_xyz_form; }
  bool use_xyz_form () const { return xyz_form_; }
};

// From Lauritzen et al, A standard test case suite for two-dimensional linear
// transport on the sphere, Geosci. Model Dev., 2012.
class InitialCondition {
  static const char* inputs[];

  static inline Real GH (const Real x, const Real y, const Real z,
                         const Real xi, const Real yi, const Real zi) {
    const Real h_max = 0.95, b = 5;
    return h_max*std::exp(-b*slmm::square(x - xi) + slmm::square(y - yi) +
                          slmm::square(z - zi));
  }

  static inline Real CB (const Real ri, const Real r) {
    const Real h_max = 1;
    return 0.5*h_max*(1 + std::cos(M_PI*ri/r));
  }

public:
  enum Shape {
    XYZTrig, GaussianHills, CosineBells, SlottedCylinders,
    CorrelatedCosineBells
  };

  static Shape from_string (const std::string& si) {
    std::string s(si);
    slmm::tolower(s);
    if (s == inputs[0]) return XYZTrig;
    if (s == inputs[1]) return GaussianHills;
    if (s == inputs[2]) return CosineBells;
    if (s == inputs[3]) return SlottedCylinders;
    if (s == inputs[4]) return CorrelatedCosineBells;
    throw std::runtime_error(si + " is not an initial condition.");
  }

  static void init (const Shape shape, const Size n, const Real* const lat,
                    const Real* const lon, Real* const u) {
    const Real lon1 = 5*(M_PI/6), lat1 = 0, lon2 = 7*(M_PI/6), lat2 = 0;
    Real x1, y1, z1, x2, y2, z2;
    slmm::ll2xyz(lat1, lon1, x1, y1, z1);
    slmm::ll2xyz(lat2, lon2, x2, y2, z2);
    switch (shape) {
    case XYZTrig: {
      for (Size i = 0; i < n; ++i) {
        Real x, y, z;
        slmm::ll2xyz(lat[i], lon[i], x, y, z, 1);
        u[i] = std::sin(3*x)*std::sin(3*y)*std::sin(4*z);
      }
    } break;
    case GaussianHills: {
      for (Size i = 0; i < n; ++i) {
        Real x, y, z;
        slmm::ll2xyz(lat[i], lon[i], x, y, z, 1);
        u[i] = GH(x, y, z, x1, y1, z1) + GH(x, y, z, x2, y2, z2);
      }
    } break;
    case CosineBells: {
      const Real r = 0.5, b = 0.1, c = 0.9;
      for (Size i = 0; i < n; ++i) {
        const Real r1 = slmm::great_circle_dist(lat[i], lon[i], lat1, lon1);
        Real h = 0;
        if (r1 < r)
          h = CB(r1, r);
        else {
          const Real r2 = slmm::great_circle_dist(lat[i], lon[i], lat2, lon2);
          if (r2 < r)
            h = CB(r2, r);
        }
        u[i] = b + c*h;
      }
    } break;
    case SlottedCylinders: {
      const Real b = 0.1, c = 1, R = 1, r = 0.5*R, lon_thr = r/(6*R),
        lat_thr = 5*(r/(12*R));
      for (Size i = 0; i < n; ++i) {
        const Real r1 = slmm::great_circle_dist(lat[i], lon[i], lat1, lon1);
        if (r1 <= r) {
          if (std::abs(lon[i] - lon1) >= lon_thr) {
            u[i] = c;
            continue;
          }
          if (std::abs(lon[i] - lon1) < lon_thr && lat[i] - lat1 < -lat_thr) {
            u[i] = c;
            continue;
          }
        }
        const Real r2 = slmm::great_circle_dist(lat[i], lon[i], lat2, lon2);
        if (r2 <= r) {
          if (std::abs(lon[i] - lon2) >= lon_thr) {
            u[i] = c;
            continue;
          }
          if (std::abs(lon[i] - lon2) < lon_thr && lat[i] - lat2 > lat_thr) {
            u[i] = c;
            continue;
          }
        }
        u[i] = b;
      }      
    } break;
    case CorrelatedCosineBells: {
      const Real a = -0.8, b = 0.9;
      init(CosineBells, n, lat, lon, u);
      for (Size i = 0; i < n; ++i)
        u[i] = a*slmm::square(u[i]) + b;
    } break;
    default: assert(0);
    }
  }

  static std::string get_inputs ()
  { return slmm::format_strings_as_list(inputs, 5); }
};

// Convert from (u,v), where u is velocity along latitude and v is velocity
// along longitude, to (x,y,z), which is velocity in the global cartesian
// coordinate system. Add a w (local vertical) component to push the position
// (X,Y,Z) back to the unit sphere.
inline void uv2xyz (
  const Real X, const Real Y, const Real Z, // position
  const Real u, const Real v, // velocity in tangent coord system
  Real& x, Real& y, Real& z) // velocity in global coord system
{
  // r should be 1 but will numerically drift, so measure it ...
  const Real r = std::sqrt(X*X + Y*Y + Z*Z);
  // ... and then add a local vertical velocity to project back to the sphere.
  const Real w = (1 - r)/slmm::consts::earth_radius_m;
  Real R[9]; // Row major.
  // The local vertical is just the position vector.
  R[2] = X/r; R[5] = Y/r; R[8] = Z/r;
  // The local along-latitude vector.
  R[0] = -Y; R[3] = X; R[6] = 0;
  const Real den = std::sqrt(R[0]*R[0] + R[3]*R[3]);
  R[0] /= den; R[3] /= den;
  // Local vertical x along-latitude.
  R[1] = R[5]*R[6] - R[8]*R[3];
  R[4] = R[8]*R[0] - R[2]*R[6];
  R[7] = R[2]*R[3] - R[5]*R[0];
  // Transform.
  x = R[0]*u + R[1]*v + R[2]*w;
  y = R[3]*u + R[4]*v + R[5]*w;
  z = R[6]*u + R[7]*v + R[8]*w;
}

// Integrate the ODE in lat-lon space. Not good numerically in the lon direction
// because of the poles.
struct Dcmip1d3llOdeFn : public OdeFnBasicRecorder {
  bool eval (const Real t, const Real* const d, Real* const f) const {
    assert ( ! use_xyz_form());
    const Real
      a = M_PI/6,
      a_ref = slmm::consts::earth_radius_m,
      tau = 1036800,
      u0 = 2*M_PI*a_ref/tau,
      sina = std::sin(a),
      cosa = std::sqrt(1 - slmm::square(sina)),
      lat = d[0],
      lon = d[1],
      sinp = std::sin(lat),
      cosp = std::cos(lat),
      sinl = std::sin(lon),
      cosl = std::cos(lon);
    // In what follows,
    //     u = u0*(cosp*cosa + sinp*cosl*sina)
    //     v = -u0*sinl*sina
    //     w = 0
    //     lat_t = slmm::m2radlat(v)
    //     lon_t = slmm::m2radlon(lat, u).
    // For numerical reasons, write this a little differently.
    const Real v = -u0*sinl*sina;
    f[0] = slmm::m2radlat(v);
    // tan(phi) is singular at the pole. We could introduce a cutoff so the wind
    // speed is not infinite, but for now it does not matter.
    f[1] = slmm::m2radlat(u0*(slmm::sign(cosp)*cosa +
                              sinp*cosl*sina/std::abs(cosp)));
    return true;
  }
};

// Also from Lauritzen et al.
struct NonDivergentWindField : public OdeFnBasicRecorder {
  bool eval (const Real t, const Real* const d, Real* const f) const {
    Real theta, lambda;
    if (use_xyz_form())
      xyz2ll(d[0], d[1], d[2], theta, lambda);
    else {
      theta = d[0];  // latitude
      lambda = d[1]; // longitude
    }
    const Real
      T = slmm::day2sec(12),
      R = slmm::consts::earth_radius_m,
      lambda_p = lambda - 2*M_PI*t/T,
      costh = std::cos(theta),
      cost = std::cos(M_PI*t/T);
    // v
    f[0] = 10*R/T*std::sin(2*lambda_p)*costh*cost;
    // u
    f[1] = R/T*(10*slmm::square(std::sin(lambda_p))*std::sin(2*theta)*cost +
                2*M_PI*costh);
    if (use_xyz_form())
      uv2xyz(d[0], d[1], d[2], f[1]/R, f[0]/R, f[0], f[1], f[2]);
    else {  
      f[0] = slmm::m2radlat(f[0]);
      f[1] = slmm::m2radlon(theta, f[1]);
    }
    return true;
  }
};

// Also from Lauritzen et al.
struct DivergentWindField : public OdeFnBasicRecorder {
  bool eval (const Real t, const Real* const d, Real* const f) const {
    Real theta, lambda;
    if (use_xyz_form())
      xyz2ll(d[0], d[1], d[2], theta, lambda);
    else {
      theta = d[0];  // latitude
      lambda = d[1]; // longitude
    }
    const Real
      T = slmm::day2sec(12),
      R = slmm::consts::earth_radius_m,
      lambda_p = lambda - 2*M_PI*t/T,
      costh = std::cos(theta),
      cost = std::cos(M_PI*t/T);
    // v
    f[0] = 2.5*R/T*std::sin(lambda_p)*slmm::cube(costh)*cost;
    // u
    f[1] = R/T*(-5*slmm::square(std::sin(0.5*lambda_p))*std::sin(2*theta)*
                slmm::square(costh)*cost + 2*M_PI*costh);
    if (use_xyz_form())
      uv2xyz(d[0], d[1], d[2], f[1]/R, f[0]/R, f[0], f[1], f[2]);
    else {  
      f[0] = slmm::m2radlat(f[0]);
      f[1] = slmm::m2radlon(theta, f[1]);
    }
    return true;
  }
};

struct NonDivergentWindFieldHack : public OdeFnBasicRecorder {
  bool eval (const Real t, const Real* const d, Real* const f) const {
    Real theta, lambda;
    if (use_xyz_form())
      xyz2ll(d[0], d[1], d[2], theta, lambda);
    else {
      theta = d[0];  // latitude
      lambda = d[1]; // longitude
    }
    const Real
      T = slmm::day2sec(12),
      R = slmm::consts::earth_radius_m,
      lambda_p = lambda,
      costh = std::cos(theta),
      cost = std::cos(M_PI*t/T);
    // v
    f[0] = 10*R/T*std::sin(2*lambda_p)*costh*cost;
    // u
    f[1] = 10*R/T*slmm::square(std::sin(lambda_p))*std::sin(2*theta)*cost;
    if (use_xyz_form())
      uv2xyz(d[0], d[1], d[2], f[1]/R, f[0]/R, f[0], f[1], f[2]);
    else {  
      f[0] = slmm::m2radlat(f[0]);
      f[1] = slmm::m2radlon(theta, f[1]);
    }
    return true;
  }
};

struct WindFieldType {
  static const char* inputs[];
public:
  enum Enum { Dcmip1d3ll, NonDivergentWindField, DivergentWindField, Rotate,
              NonDivergentWindFieldHack };
  static Enum from_string (const std::string& si) {
    std::string s(si);
    slmm::tolower(s);
    if (s == inputs[0]) return Dcmip1d3ll;
    if (s == inputs[1]) return NonDivergentWindField;
    if (s == inputs[2]) return DivergentWindField;
    if (s == inputs[3]) return Rotate;
    if (s == inputs[4]) return NonDivergentWindFieldHack;
    throw std::runtime_error(si + " is not an ODE function.");
  }
  static std::string get_inputs ()
  { return slmm::format_strings_as_list(inputs, 4); }
};

} // namespace gallery
} // namespace slmm

#endif
