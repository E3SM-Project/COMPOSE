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
    SmoothBelts, CBandSC, Zero, EquatorStep, EquatorSmoothStep
  };

  static Shape from_string(const std::string& si);
  static const char* to_string(const Shape& shape);

  // Input is in radians.
  static void init(const Shape shape, const Size n, const Real* const lat,
                   const Real* const lon, Real* const u);

  static std::string get_inputs();
};

// Integrate the ODE in lat-lon space. Not good numerically in the lon direction
// because of the poles.
struct Dcmip1d3llOdeFn : public OdeFn {
  bool eval(const Real t, const Real* const d, Real* const f) const override;
};

// Also from Lauritzen et al.
struct NonDivergentWindField : public OdeFn {
  bool eval(const Real t, const Real* const d, Real* const f) const override;
};

// Also from Lauritzen et al.
struct DivergentWindField : public OdeFn {
  bool eval(const Real t, const Real* const d, Real* const f) const override;
};

struct MovingVortices : public OdeFn {
  bool eval(const Real t, const Real* const d, Real* const f) const override;
  static const Real rho0, gamma;
  // Omega is the solid-body rotation rate.
  static Real get_Omega();
  static Real calc_rho(const Real theta, const Real lambda);
  // omega is the vortex rotation.
  static Real calc_omega(const Real Omega, const Real rho);
  // Time is in seconds.
  static void calc_tracer(const Real time, const Size n, const Real* const lat,
                          const Real* const lon, Real* const u);
};

struct NonDivergentWindFieldHack : public OdeFn {
  bool eval(const Real t, const Real* const d, Real* const f) const override;
};

// DCMIP test 1-1 with w=0 and a particular value of z.
struct TestWindField : public OdeFn {
  bool eval(const Real t, const Real* const d, Real* const f) const override;
};

struct WindFieldType {
  static const char* inputs[];
public:
  enum Enum { Dcmip1d3ll, NonDivergentWindField, DivergentWindField, Rotate,
              NonDivergentWindFieldHack, MovingVortices, TestWindField };
  static Enum from_string (const std::string& si) {
    std::string s(si);
    slmm::tolower(s);
    if (s == inputs[0]) return Dcmip1d3ll;
    if (s == inputs[1]) return NonDivergentWindField;
    if (s == inputs[2]) return DivergentWindField;
    if (s == inputs[3]) return Rotate;
    if (s == inputs[4]) return NonDivergentWindFieldHack;
    if (s == inputs[5]) return MovingVortices;
    if (s == inputs[6]) return TestWindField;
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
