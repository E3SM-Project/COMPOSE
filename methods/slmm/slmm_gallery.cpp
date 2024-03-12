#include "slmm_gallery.hpp"

namespace slmm {
namespace gallery {

const char* InitialCondition::inputs[] =
  {"xyztrig", "gaussianhills", "cosinebells", "slottedcylinders",
   "correlatedcosinebells", "constant", "vortextracer", "toychem1", "toychem2",
   "slotcyltrig", "smoothbelts", "cbandsc", "zero", "equatorstep",
   "equatorsmoothstep"};

const char* WindFieldType::inputs[] =
  {"dcmip1d3ll", "nondivergent", "divergent", "rotate", "nondivergenthack",
   "movingvortices", "testfn"};

std::string InitialCondition::get_inputs () {
  return slmm::format_strings_as_list(inputs, 10);
}

InitialCondition::Shape InitialCondition::from_string (const std::string& si) {
  std::string s(si);
  slmm::tolower(s);
  if (s == inputs[ 0]) return XYZTrig;
  if (s == inputs[ 1]) return GaussianHills;
  if (s == inputs[ 2]) return CosineBells;
  if (s == inputs[ 3]) return SlottedCylinders;
  if (s == inputs[ 4]) return CorrelatedCosineBells;
  if (s == inputs[ 5]) return Constant;
  if (s == inputs[ 6]) return VortexTracer;
  if (s == inputs[ 7]) return ToyChem1;
  if (s == inputs[ 8]) return ToyChem2;
  if (s == inputs[ 9]) return SlotCylTrig;
  if (s == inputs[10]) return SmoothBelts;
  if (s == inputs[11]) return CBandSC;
  if (s == inputs[12]) return Zero;
  if (s == inputs[13]) return EquatorStep;
  if (s == inputs[14]) return EquatorSmoothStep;
  throw std::runtime_error(si + " is not an initial condition.");
}

const char* InitialCondition::to_string (const Shape& shape) {
  switch (shape) {
  case XYZTrig: return inputs[0];
  case GaussianHills: return inputs[1];
  case CosineBells: return inputs[2];
  case SlottedCylinders: return inputs[3];
  case CorrelatedCosineBells: return inputs[4];
  case Constant: return inputs[5];
  case VortexTracer: return inputs[6];
  case ToyChem1: return inputs[7];
  case ToyChem2: return inputs[8];
  case SlotCylTrig: return inputs[9];
  case SmoothBelts: return inputs[10];
  case CBandSC: return inputs[11];
  case Zero: return inputs[12];
  case EquatorStep: return inputs[13];
  case EquatorSmoothStep: return inputs[14];
  }
  throw std::runtime_error("Should not be here.");
}

static inline Real GH (const Real x, const Real y, const Real z,
                       const Real xi, const Real yi, const Real zi) {
  const Real h_max = 0.95, b = 5;
  const Real r2 = (slmm::square(x - xi) + slmm::square(y - yi) +
                   slmm::square(z - zi));
  return h_max*std::exp(-b*r2);
}

static inline Real CB (const Real ri, const Real r) {
  const Real h_max = 1;
  return 0.5*h_max*(1 + std::cos(M_PI*ri/r));
}

void InitialCondition::init (const Shape shape, const Size n, const Real* const lat,
                             const Real* const lon, Real* const u) {
  const Real lon1 = 5*(M_PI/6), lat1 = 0, lon2 = -5*(M_PI/6), lat2 = 0;
  Real x1, y1, z1, x2, y2, z2;
  slmm::ll2xyz(lat1, lon1, x1, y1, z1);
  slmm::ll2xyz(lat2, lon2, x2, y2, z2);
  switch (shape) {
  case XYZTrig: {
    for (Size i = 0; i < n; ++i) {
      Real x, y, z;
      slmm::ll2xyz(lat[i], lon[i], x, y, z, 1);
      u[i] = 0.5*(1 + std::sin(3*x)*std::sin(3*y)*std::sin(4*z));
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
  case VortexTracer: {
    MovingVortices::calc_tracer(0, n, lat, lon, u);
  } break;
  case CorrelatedCosineBells: {
    const Real a = -0.8, b = 0.9;
    init(CosineBells, n, lat, lon, u);
    for (Size i = 0; i < n; ++i)
      u[i] = a*slmm::square(u[i]) + b;
  } break;
  case Constant: {
    for (Size i = 0; i < n; ++i)
      u[i] = 0.42;
  } break;
  case Zero: {
    for (Size i = 0; i < n; ++i)
      u[i] = 0;
  } break;
  case EquatorStep: {
    for (Size i = 0; i < n; ++i)
      u[i] = lat[i] >= 0 ? 1 : 0.1;
  } break;
  case EquatorSmoothStep: {
    static const Real lat_thr = M_PI/4, a = 0.1, b = 1;
    for (Size i = 0; i < n; ++i) {
      if (std::abs(lat[i]) >= lat_thr)
        u[i] = lat[i] >= 0 ? b : a;
      else
        u[i] = a + ((b-a)/2)*(1 + std::sin(M_PI/2*(lat[i]/lat_thr)));
    }
  } break;
  case ToyChem1: {
    for (Size i = 0; i < n; ++i) {
      Real k1, k2;
      ToyChem::k_vals(lat[i], lon[i], k1, k2);
      const Real r = k1 / (4*k2);
      const Real det = std::sqrt(r*r + 2*ToyChem::constant*r);
      u[i] = det - r;
    }
  } break;
  case ToyChem2: {
    init(ToyChem1, n, lat, lon, u);
    for (Size i = 0; i < n; ++i)
      u[i] = ToyChem::constant - u[i];
  } break;
  case SlotCylTrig: {
    init(SlottedCylinders, n, lat, lon, u);
    for (Size i = 0; i < n; ++i) {
      if (u[i] > 0.9) continue;
      Real x, y, z;
      slmm::ll2xyz(lat[i], lon[i], x, y, z, 1);
      u[i] = 0.5*(1 + std::sin(3*x)*std::sin(3*y)*std::sin(4*z));
    }
  } break;
  case SmoothBelts: {
    static bool first = true;
    static Real R[9];
    if (first) {
      const Real axis[] = {1, 0, 0};
      const Real angle = 0.1*M_PI;
      form_rotation(axis, angle, R);
      first = false;
    }
    for (Size i = 0; i < n; ++i) {
      Real x[3];
      slmm::ll2xyz(lat[i], lon[i], x[0], x[1], x[2], 1);
      const Real y2 = R[6]*x[0] + R[7]*x[1] + R[8]*x[2];
      u[i] = 0.5*(1 + std::cos(M_PI*y2));
    }
  } break;
  case CBandSC: {
    for (Size i = 0; i < n; ++i) {
      {
        const Real lon1 = 3*(M_PI/6);
        const Real r = 0.5, b = 0.1, c = 0.9;
        u[i] = b;
        const Real r1 = slmm::great_circle_dist(lat[i], lon[i], lat1, lon1);
        if (r1 < r) u[i] += c*CB(r1, r);
      }
      {
        const Real lon2 = -5*(M_PI/6);
        const Real c = 0.9, R = 1, r = 0.5*R, lon_thr = r/(6*R),
          lat_thr = 5*(r/(12*R));
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
      }
    }
  } break;
  default: assert(0);
  }
}

void ToyChem::k_vals (const Real& lat, const Real& lon, Real& k1, Real& k2) {
  k1 = std::max<Real>(0, (std::sin(lat)*std::sin(k1_lat_center) +
                          (std::cos(lat)*std::cos(k1_lat_center)*
                           std::cos(lon - k1_lon_center))));
  k2 = 1;
}

void ToyChem::tendency (const Real& lat, const Real& lon,
                        Real cl, Real cl2, // mixing ratio mass/mass
                        const Real& dt, Real& cl_f, Real& cl2_f) {
  cl2 = 0.5*cl2; // convert to atom number

  Real k1, k2;
  k_vals(lat, lon, k1, k2);
  const Real r = k1 / (4*k2);
  const Real cly = cl + 2*cl2;

  const Real det = std::sqrt(r*r + 2*r*cly);
  const Real expdt = std::exp(-4*k2*det*dt);

  Real el;
  if (std::abs(det * k2 * dt) > 1e-16)
    el = (1 - expdt)/(det*dt);
  else
    el = 4*k2;

  cl_f  = -el * (cl - det + r)*(cl + det + r) / (1 + expdt + dt*el*(cl + r));
  cl2_f = -cl_f; // no /2 b/c there's an implicit conversion back to mixing ratio
}

// Convert from (u,v), where u is velocity along latitude and v is velocity
// along longitude, to (x,y,z), which is velocity in the global cartesian
// coordinate system. Add a w (local vertical) component to push the position
// (X,Y,Z) back to the unit sphere.
static void uv2xyz (
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

bool Dcmip1d3llOdeFn
::eval (const Real t, const Real* const d, Real* const f) const {
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

bool NonDivergentWindField
::eval (const Real t, const Real* const d, Real* const f) const {
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

bool DivergentWindField
::eval (const Real t, const Real* const d, Real* const f) const {
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

// This test is based on
//     Nair, R. D., & Jablonowski, C. (2008). Moving vortices on the sphere: A
//     test case for horizontal advection problems. Monthly Weather Review,
//     136(2), 699-711.
// The code uses the formulas in
//     Bosler, P. A., Kent, J., Krasny, R., & Jablonowski, C. (2017). A
//     Lagrangian particle method with remeshing for tracer transport on the
//     sphere. Journal of Computational Physics, 340, 639-654.

const Real MovingVortices::rho0 = 3;
const Real MovingVortices::gamma = 5;

Real MovingVortices::get_Omega () {
  return 2*M_PI/slmm::day2sec(12);
}

Real MovingVortices::calc_rho (const Real theta, const Real lambda) {
  return rho0*std::sqrt(1 - slmm::square(std::cos(theta)*std::sin(lambda)));
}

Real MovingVortices::calc_omega (const Real Omega, const Real rho) {
  return (rho == 0 ?
          0 :
          (Omega*slmm::consts::earth_radius_m*
           1.5*std::sqrt(3)*std::tanh(rho) /
           (rho*slmm::square(std::cosh(rho)))));
}

void MovingVortices
::calc_tracer (const Real time, const Size n, const Real* const lat,
               const Real* const lon, Real* const u) {
  const auto Omega = get_Omega(), R = slmm::consts::earth_radius_m;
  for (Size i = 0; i < n; ++i) {
    const Real
      lon_d = lon[i] - Omega*time,
      lambda_p = std::atan2(-std::cos(lon_d), std::tan(lat[i])),
      rho = calc_rho(lat[i], lon_d),
      omega = calc_omega(Omega, rho);
    u[i] = 1 - std::tanh((rho/MovingVortices::gamma)*
                         std::sin(lambda_p - (omega/R)*time));
  }
}

bool MovingVortices
::eval (const Real t, const Real* const d, Real* const f) const {
  Real theta, lambda;
  if (use_xyz_form())
    xyz2ll(d[0], d[1], d[2], theta, lambda);
  else {
    theta = d[0];
    lambda = d[1];
  }
  const Real
    R = slmm::consts::earth_radius_m,
    Omega = get_Omega(),
    // Solid body rotation factor. 1 is one rotation per T. 0 gives stationary
    // vortices. 1e-3 is enough to stabilize the instability of the stationary
    // vortex center being on a cell corner.
    fac = 1,
    lambda_p = lambda - fac*Omega*t,
    costh = std::cos(theta),
    rho = calc_rho(theta, lambda_p),
    omega = calc_omega(Omega, rho);
  // v
  f[0] = omega*std::cos(lambda_p);
  // u
  f[1] = omega*std::sin(lambda_p)*std::sin(theta) + R*Omega*costh*fac;
  if (use_xyz_form())
    uv2xyz(d[0], d[1], d[2], f[1]/R, f[0]/R, f[0], f[1], f[2]);
  else {
    f[0] = slmm::m2radlat(f[0]);
    f[1] = slmm::m2radlon(theta, f[1]);
  }
  return true;
}

bool NonDivergentWindFieldHack
::eval (const Real t, const Real* const d, Real* const f) const {
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

bool TestWindField
::eval (const Real t, const Real* const d, Real* const f) const {
  Real theta, lambda;
  if (use_xyz_form())
    xyz2ll(d[0], d[1], d[2], theta, lambda);
  else {
    theta = d[0];  // latitude
    lambda = d[1]; // longitude
  }
  const Real
    ztop = 12000,
    z = 0.05*ztop,
    T = slmm::day2sec(12),
    R = slmm::consts::earth_radius_m,
    T0 = 300,
    Rd = 287.04,
    g  = 9.80616,
    p0 = 100000,
    H  = Rd*T0/g,
    u0 = 2*M_PI*R/T,
    k0 = 10*R/T,
    omega0 = (2*23000*M_PI)/T,
    lambda_p = lambda - 2*M_PI*t/T,
    costh   = std::cos(theta),
    cost    = std::cos(M_PI*t/T),
    sinlamp = std::sin(lambda_p),
    p    = p0*std::exp(-z/H),
    ptop = p0*std::exp(-ztop/H),
    bs = 0.2;

  const Real
    s = (1.0 + std::exp((ptop-p0)/(bs*ptop)) - std::exp((p-p0)/(bs*ptop)) -
         std::exp((ptop-p)/(bs*ptop))),
    s_p = (-std::exp((p-p0)/(bs*ptop)) + std::exp((ptop-p)/(bs*ptop)))/(bs*ptop),
    ud = (omega0*R)*std::cos(lambda_p)*square(costh)*cost*s_p;

  // v
  f[0] = 10*R/T*std::sin(2*lambda_p)*costh*cost;
  // u
  f[1] = R/T*(10*slmm::square(std::sin(lambda_p))*std::sin(2*theta)*cost +
              2*M_PI*costh) + ud;

  if (use_xyz_form())
    uv2xyz(d[0], d[1], d[2], f[1]/R, f[0]/R, f[0], f[1], f[2]);
  else {  
    f[0] = slmm::m2radlat(f[0]);
    f[1] = slmm::m2radlon(theta, f[1]);
  }
  return true;
}

} // namespace gallery
} // namespace slmm
