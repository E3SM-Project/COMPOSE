#include "slmm_gallery.hpp"

namespace slmm {
namespace gallery {

const char* InitialCondition::inputs[] =
  {"xyztrig", "gaussianhills", "cosinebells", "slottedcylinders",
   "correlatedcosinebells", "constant", "vortextracer", "toychem1", "toychem2",
   "slotcyltrig", "smoothbelts", "cbandsc", "zero"};

const char* WindFieldType::inputs[] =
  {"dcmip1d3ll", "nondivergent", "divergent", "rotate", "nondivergenthack",
   "movingvortices"};

std::string InitialCondition::get_inputs () {
  return slmm::format_strings_as_list(inputs, 10);
}

InitialCondition::Shape InitialCondition::from_string (const std::string& si) {
  std::string s(si);
  slmm::tolower(s);
  if (s == inputs[0]) return XYZTrig;
  if (s == inputs[1]) return GaussianHills;
  if (s == inputs[2]) return CosineBells;
  if (s == inputs[3]) return SlottedCylinders;
  if (s == inputs[4]) return CorrelatedCosineBells;
  if (s == inputs[5]) return Constant;
  if (s == inputs[6]) return VortexTracer;
  if (s == inputs[7]) return ToyChem1;
  if (s == inputs[8]) return ToyChem2;
  if (s == inputs[9]) return SlotCylTrig;
  if (s == inputs[10]) return SmoothBelts;
  if (s == inputs[11]) return CBandSC;
  if (s == inputs[12]) return Zero;
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
    for (Size i = 0; i < n; ++i) {
      const Real lambda_p = std::atan2(-std::cos(lon[i]), std::tan(lat[i]));
      const Real rr = 3*std::sqrt((1 - slmm::square(std::cos(lat[i])) *
                                   slmm::square(std::sin(lon[i]))));
      u[i] = 0.5*(1 - std::tanh(0.2*rr*std::sin(lambda_p)));
    }
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

} // namespace gallery
} // namespace slmm
