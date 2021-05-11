#ifndef INCLUDE_SLMM_GALLERY_HPP
#define INCLUDE_SLMM_GALLERY_HPP

#include "slmm_defs.hpp"
#include "slmm_time_int.hpp"

namespace slmm {
namespace gallery {

class OdeFn {
  mutable int ne_;
  bool xyz_form_;
public:
  OdeFn () : ne_(0), xyz_form_(false) {}
  virtual bool eval (const Real t, const Real* const d, Real* const f) const = 0;
  virtual void record (const Real t, const Real* const y) const { ++ne_; }
  int ne () const { return ne_; }
  void set_xyz_form (const bool use_xyz_form) { xyz_form_ = use_xyz_form; }
  bool use_xyz_form () const { return xyz_form_; }
};

// From Lauritzen et al, A standard test case suite for two-dimensional linear
// transport on the sphere, Geosci. Model Dev., 2012.
class InitialCondition {
  static const char* inputs[];

public:
  enum Shape {
    XYZTrig, GaussianHills, CosineBells, SlottedCylinders,
    CorrelatedCosineBells, Constant, VortexTracer, ToyChem1, ToyChem2,
    SlotCylTrig, // Slotted cylinders with a trig background
    SmoothBelts, CBandSC, Zero
  };

  static Shape from_string(const std::string& si);
  static const char* to_string(const Shape& shape);

  // Input is in radians.
  static void init(const Shape shape, const Size n, const Real* const lat,
                   const Real* const lon, Real* const u);

  static std::string get_inputs();
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
struct Dcmip1d3llOdeFn : public OdeFn {
  bool eval (const Real t, const Real* const d, Real* const f) const override {
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
struct NonDivergentWindField : public OdeFn {
  bool eval (const Real t, const Real* const d, Real* const f) const override {
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
struct DivergentWindField : public OdeFn {
  bool eval (const Real t, const Real* const d, Real* const f) const override {
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

struct MovingVortices : public OdeFn {
  bool eval( const Real t, const Real* const d, Real* const f) const override {
    Real theta, lambda;
    if (use_xyz_form()) {
      xyz2ll(d[0], d[1], d[2], theta, lambda);
    }
    else {
      theta = d[0];
      lambda = d[1];
    }
    const Real T = slmm::day2sec(12);           // rotational period, in seconds
    const Real R = slmm::consts::earth_radius_m;// sphere radius, in meters
    const Real Omega = 2*M_PI/T;                // angular velocity of background rotation, radians per second
    const Real u0 = Omega*R;                    // velocity due to background rotation, meters per second
    const Real lambda_p = lambda - Omega*t;     // longitudinal center of vortex, radians
    const Real costh = std::cos(theta);         // cosine of latitude
        
    const Real rr = 3*std::sqrt(1 - slmm::square(costh)*slmm::square(std::sin(lambda_p)));
    const Real rr_denom = rr / (slmm::square(rr) + slmm::square(1.0e-10));
    const Real omg = 1.5*std::sqrt(3.0)*u0*std::tanh(rr) * rr_denom / (slmm::square(std::cosh(rr)));
    // v
    f[0] = omg * std::cos(lambda_p);
    // u
    f[1] = omg * std::sin(lambda_p) * std::sin(theta) + u0 * costh;
    if (use_xyz_form())
      uv2xyz(d[0], d[1], d[2], f[1]/R, f[0]/R, f[0], f[1], f[2]);
    else {
      f[0] = slmm::m2radlat(f[0]);
      f[1] = slmm::m2radlon(theta, f[1]);
    }
    return true;
  }
};

struct NonDivergentWindFieldHack : public OdeFn {
  bool eval (const Real t, const Real* const d, Real* const f) const override {
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
              NonDivergentWindFieldHack, MovingVortices };
  static Enum from_string (const std::string& si) {
    std::string s(si);
    slmm::tolower(s);
    if (s == inputs[0]) return Dcmip1d3ll;
    if (s == inputs[1]) return NonDivergentWindField;
    if (s == inputs[2]) return DivergentWindField;
    if (s == inputs[3]) return Rotate;
    if (s == inputs[4]) return NonDivergentWindFieldHack;
    if (s == inputs[5]) return MovingVortices;
    throw std::runtime_error(si + " is not an ODE function.");
  }
  static std::string get_inputs ()
  { return slmm::format_strings_as_list(inputs, 6); }
};

// From toy chemistry model ref impl chemistry.F90.
struct ToyChem {
  static constexpr Real deg2rad = M_PI/180;
  static constexpr Real constant = 4e-6;
  static constexpr Real k1_lat_center =  20*deg2rad;
  static constexpr Real k1_lon_center = 300*deg2rad;

  // Input is in radians.
  static void k_vals(const Real& lat, const Real& lon, Real& k1, Real& k2);

  // Input is in radians.
  static void tendency(const Real& lat, const Real& lon,
                       Real cl, Real cl2, // mixing ratio mass/mass
                       const Real& dt, Real& cl_f, Real& cl2_f);
};

} // namespace gallery
} // namespace slmm

#endif
