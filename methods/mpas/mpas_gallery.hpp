#ifndef INCLUDE_MPAS_GALLERY_HPP
#define INCLUDE_MPAS_GALLERY_HPP

#include "mpas_mesh.hpp"
#include "mpas_comm.hpp"

// From slmm::gallery.

// From Lauritzen et al, A standard test case suite for two-dimensional linear
// transport on the sphere, Geosci. Model Dev., 2012.
class InitialCondition {
  static const char* inputs[];

public:
  enum Shape {
    XYZTrig, GaussianHills, CosineBells, SlottedCylinders,
    CorrelatedCosineBells, Constant, VortexTracer, ToyChem1, ToyChem2,
    SlotCylTrig // Slotted cylinders with a trig background
  };

  static Shape from_string(const std::string& si);
  static const char* to_string(const Shape& shape);

  // Input is in radians.
  static void init(const Shape shape, const Size n, const Real* const lon,
                   const Real* const lat, Real* const u);

  static std::string get_inputs();
};

// Get the IC at cell ci.
Real get_ic(const MpasMesh& m, const InitialCondition::Shape ic, const int ci);

Real get_ic_dcmip12t11(const MpasMesh& m, const int qi, const int ci,
                       const Real zeta0, const Real zeta1);

struct CmpOutput { Real mass, l1, l2, linf, lo, hi; };

void cmp(const MpasMesh& m, Comm& he, const Int qi, const Int lvl, const Int nlvl,
         const Int ncown,
         // if rho is provided, use q rather than Q in error norms
         const Ar2* rhotrue, const Ar3& Qtrue, const Ar2* rho, const Ar3& Q,
         // optionally provide cell area
         CmpOutput& out, const Real* cell_area = nullptr);

struct VelocityField {
  typedef std::shared_ptr<VelocityField> Ptr;

  virtual ~VelocityField () {}

  virtual void operator() (const Real t, const Real lon, const Real lat,
                           const Real z, const Int level, // z and level may be unused
                           Real& u, Real& v, Real& w) const = 0;
  virtual std::string name () const { return ""; }
  virtual Real cycle_time () const { return 12*86400; }
  // Loose upper bound is fine.
  virtual Real max_vel() const;
  virtual bool is_3d () const { return false; }
  // Interface level (ilvl) and midpoint level (mlvl) to z coordinate, rho,
  // shape profile.
  virtual Real get_ic_shape_fn(const Real z0, const Real z1) const;
  virtual Real get_rho (const Real z0, const Real z1) const { return 1; }
  // Default levels. But a VelocityField impl should not rely on these being
  // used except for the top and bottom values.
  virtual Real mlvl_to_z(const int nmlvl, const int mlvl) const;
  virtual Real ilvl_to_z (const int nilvl, const int ilvl) const { return 0; }
};

struct Rotation : public VelocityField {
  Rotation(bool zaxis = false);
  void operator()(const Real t, const Real lon, const Real lat, const Real z,
                  const Int level, Real& u, Real& v, Real& w) const override;
  std::string name () const override { return "rotation"; }
private:
  Vec Omega;
};

struct NonDivergentWindField : public VelocityField {
  void operator()(const Real t, const Real lon, const Real lat, const Real,
                  const Int, Real& u, Real& v, Real& w) const override;
  std::string name () const override { return "nondiv"; }
};

struct DivergentWindField : public VelocityField {
  void operator()(const Real t, const Real lon, const Real lat, const Real,
                  const Int, Real& u, Real& v, Real& w) const override;
  std::string name () const override { return "div"; }
};

struct Multiple2DWindFields : public VelocityField {
  void operator()(const Real t, const Real lon, const Real lat, const Real,
                  const Int level, Real& u, Real& v, Real& w) const override;
  std::string name () const override { return "mult"; }
};

// DCMIP 2012 test 1-1.
struct Dcmip20212Test11WindField : public VelocityField {
  void operator()(const Real t, const Real lon, const Real lat, const Real z,
                  const Int, Real& u, Real& v, Real& w) const override;
  std::string name () const override { return "dcmip2012t11"; }
  bool is_3d () const override { return true; }
  Real get_ic_shape_fn(const Real z0, const Real z1) const override;
  Real ilvl_to_z(const int nilvl, const int ilvl) const override;
  Real get_rho(const Real z0, const Real z1) const override;
};

// DCMIP 2012 test 1-1 with some mods.
struct NonDivergent3DWindField : public VelocityField {
  void operator()(const Real t, const Real lon, const Real lat, const Real z,
                  const Int, Real& u, Real& v, Real& w) const override;
  std::string name () const override { return "nondiv3d"; }
  bool is_3d () const override { return true; }
  Real ilvl_to_z(const int nilvl, const int ilvl) const override;
  Real get_rho(const Real z0, const Real z1) const override;
};

struct NonDivergent3DVertPWindField : public VelocityField {
  NonDivergent3DVertPWindField();
  void operator()(const Real t, const Real lon, const Real lat, const Real z,
                  const Int, Real& u, Real& v, Real& w) const override;
  std::string name () const override { return "nondiv3dvertp"; }
  bool is_3d () const override { return true; }
  Real ilvl_to_z(const int nilvl, const int ilvl) const override;
  Real get_rho(const Real z0, const Real z1) const override;
  Real get_ic_shape_fn(const Real z0, const Real z1) const override;
};

// Study DCMIP 2012 test 1-1 in one near-bdy level.
struct TestWindField : public VelocityField {
  void operator()(const Real t, const Real lon, const Real lat, const Real z,
                  const Int, Real& u, Real& v, Real& w) const override;
  Real ilvl_to_z(const int nlvl, const int lvl) const override;
  std::string name () const override { return "test"; }
};

// DCMIP 2012 test 1-2: Hadley-like meridional circulation.
struct HadleyWindField : public VelocityField {
  void operator()(const Real t, const Real lon, const Real lat, const Real z,
                  const Int, Real& u, Real& v, Real& w) const override;
  std::string name () const override { return "hadley"; }
  Real cycle_time () const override { return 86400; }
  bool is_3d () const override { return true; }
  Real get_ic_shape_fn(const Real z0, const Real z1) const override;
  Real ilvl_to_z(const int nilvl, const int ilvl) const override;
  Real get_rho(const Real z0, const Real z1) const override;
};

VelocityField::Ptr vf_from_short(char sym, const bool verbose=true);
InitialCondition::Shape ic_from_short(char sym, const bool verbose=true);

// A reasonable nstep for one cycle of the flow.
Int calc_default_nstep(const MpasMesh& m, Comm& he, const VelocityField& vf);

// Compute or adjust nhalo, nsubstep, or both based on the flow given by vf,
// with nstep the number of time steps for one cycle. Indicate parameters to set
// by assigning them < 0 on input. If neither is < 0 on input, the adjustment
// prioritizes retaining the input value of nsubstep.
void calc_substep_data(const MpasMesh& m, Comm& he, const VelocityField& vf,
                       const int nstep, int& traj_nhalo, int& nsubstep);

void set_interface_z_all_cols(const VelocityField& vf, MpasMesh& mm);

/* Define the monotonically increasing fn f(x) for x in [0,1]. Layer i in
   [0,n-1] has size
       int(df(x)/dx; i/n, (i+1)/n) = f((i+1)/n) - f(i/n);
   e.g., if f(x) = x, then layer i has size 1/n. The quotient largest over
   smallest layer is asymptotically
       max df/dx / min df/dx.
   Doubling the number of layers refines the layers by approximately 2. If f(x)
   is smooth, then so is the distortion. Finally, f(1/2) (say) is always the
   same, so as the grid is refined, the middle stays approximately fixed and
   asymptotes to f(1/2).
     We want there to be three regions, either (fine, coarse, fine) or (coarse,
   fine, coarse). This function provides the first. The center of the coarse
   region, ctr, should vary smoothly in space.
     A good f is given by the following shape fn for the domain x in [-1,1]:
       h(x) = 1/(1 + exp(-a x))
       g(x; c) = (g(c+x) - g(c-1))/(g(c+1) - g(c-1))
       f(x; ctr) = g(x; -ctr)
   ctr is where f is steepest. f(-1) = 0 and f(1) = 1. A reasonable value of a
   is 1.25, and a reasonable range for ctr is [-0.75,0.75].
     To use f, make a smooth function ctr(x,y,z), where (x,y,z) is position on
   the unit sphere. Then set zeta interfaces at f(linspace(-1,1,nilvl);
   ctr(x,y,z)).
 */
Real get_level_shape_fn(const Real x,   // in [-1,1]
                        const Real ctr, // in [-1,1] but practically in [-0.75,0.75]
                        const Real a = 1.25);

#endif
