#include "mpas_gallery.hpp"
#include "mpas_geometry.hpp"
#include "mpas_triquad.hpp"
#include "mpas_comm_util.hpp"
#include "mpas_accum.hpp"

#include <sstream>

const char* InitialCondition::inputs[] =
  {"xyztrig", "gaussianhills", "cosinebells", "slottedcylinders",
   "correlatedcosinebells", "constant", "vortextracer", "toychem1", "toychem2",
   "slotcyltrig"};

std::string format_strings_as_list (const char** strings, const Size n) {
  std::stringstream ss;
  ss << "{";
  for (Size i = 0; i < n-1; ++i) ss << strings[i] << ", ";
  ss << strings[n-1] << "}";
  return ss.str();
}

std::string format_strings_as_list(const char** strings, const Size n);

std::string& tolower (std::string& s) {
  for (auto& c: s)
    c = std::tolower(c);
  return s;
}

std::string InitialCondition::get_inputs () {
  return format_strings_as_list(inputs, 10);
}

InitialCondition::Shape InitialCondition::from_string (const std::string& si) {
  std::string s(si);
  tolower(s);
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
  }
  throw std::runtime_error("Should not be here.");
}

static inline Real GH (const Real x, const Real y, const Real z,
                       const Real xi, const Real yi, const Real zi) {
  const Real h_max = 0.95, b = 5;
  const Real r2 = (square(x - xi) + square(y - yi) +
                   square(z - zi));
  return h_max*std::exp(-b*r2);
}

static inline Real CB (const Real ri, const Real r) {
  const Real h_max = 1;
  return 0.5*h_max*(1 + std::cos(M_PI*ri/r));
}

void InitialCondition::init (const Shape shape, const Size n, const Real* const lon,
                             const Real* const lat, Real* const u) {
  const Real lon1 = 5*(M_PI/6), lat1 = 0, lon2 = -5*(M_PI/6), lat2 = 0;
  Real x1, y1, z1, x2, y2, z2;
  lonlat2xyz(lon1, lat1, x1, y1, z1);
  lonlat2xyz(lon2, lat2, x2, y2, z2);
  switch (shape) {
  case XYZTrig: {
    for (Size i = 0; i < n; ++i) {
      Real x, y, z;
      lonlat2xyz(lon[i], lat[i], x, y, z, 1);
      u[i] = 0.5*(1 + std::sin(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z));
    }
  } break;
  case GaussianHills: {
    for (Size i = 0; i < n; ++i) {
      Real x, y, z;
      lonlat2xyz(lon[i], lat[i], x, y, z, 1);
      u[i] = GH(x, y, z, x1, y1, z1) + GH(x, y, z, x2, y2, z2);
    }
  } break;
  case CosineBells: {
    const Real r = 0.5, b = 0.1, c = 0.9;
    for (Size i = 0; i < n; ++i) {
      const Real r1 = great_circle_dist(lon[i], lat[i], lon1, lat1);
      Real h = 0;
      if (r1 < r)
        h = CB(r1, r);
      else {
        const Real r2 = great_circle_dist(lon[i], lat[i], lon2, lat2);
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
      const Real r1 = great_circle_dist(lon[i], lat[i], lon1, lat1);
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
      const Real r2 = great_circle_dist(lon[i], lat[i], lon2, lat2);
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
      const Real rr = 3*std::sqrt((1 - square(std::cos(lat[i])) *
                                   square(std::sin(lon[i]))));
      u[i] = 0.5*(1 - std::tanh(0.2*rr*std::sin(lambda_p)));
    }
  } break;
  case CorrelatedCosineBells: {
    const Real a = -0.8, b = 0.9;
    init(CosineBells, n, lat, lon, u);
    for (Size i = 0; i < n; ++i)
      u[i] = a*square(u[i]) + b;
  } break;
  case Constant: {
    for (Size i = 0; i < n; ++i)
      u[i] = 0.42;
  } break;
  case SlotCylTrig: {
    init(SlottedCylinders, n, lat, lon, u);
    for (Size i = 0; i < n; ++i) {
      if (u[i] > 0.9) continue;
      Real x, y, z;
      lonlat2xyz(lon[i], lat[i], x, y, z, 1);
      u[i] = 0.5*(1 + std::sin(3*x)*std::sin(3*y)*std::sin(4*z));
    }
  } break;
  default: assert(0);
  }
}

Real get_ic (const MpasMesh& m, const InitialCondition::Shape ic, const int ci) {
#if 1
  // Inject the value at the cell center. Use this when comparing with MPAS-O
  // test cases.
  Real lon, lat;
  xyz2lonlat(m.c.x[ci], m.c.y[ci], m.c.z[ci], lon, lat);
  Real q;
  InitialCondition::init(ic, 1, &lon, &lat, &q);
  return q;
#else
  // Compute the cell mean value.
  TriangleQuadrature tq;
  const Real* tq_alpha, * tq_w;
  int tqn;
  tq.get_coef(4, tq_alpha, tq_w, tqn);
  const int ne = m.t.c2ne[ci];
  Vec v0; fill(m.c, ci, v0, true);
  Real area = 0, mass = 0;
  for (int i = 0; i < ne; ++i) {
    const int ie = m.t.get_c2e(ci,i), vi1 = m.t.get_e2v(ie,0), vi2 = m.t.get_e2v(ie,1);
    Vec v1; fill(m.v, vi1, v1, true);
    Vec v2; fill(m.v, vi2, v2, true);
    for (int qi = 0; qi < tqn; ++qi) {
      Real u[3];
      const Real jac = geo::calc_tri_jacobian(v0, v1, v2, slice(tq_alpha, qi), u);
      Real lon, lat, y;
      xyz2lonlat(u[0], u[1], u[2], lon, lat);
      InitialCondition::init(ic, 1, &lon, &lat, &y);
      area += tq_w[qi]*jac;
      mass += tq_w[qi]*jac*y;
    }
  }
  return mass/area;
#endif
}

Real get_ic_dcmip12t11 (const MpasMesh& m, const int qi, const int ci,
                        const Real zeta0, const Real zeta1) {
  assert(qi >= 0 && qi < 5);
 
  Real q = -1;
  
  if (qi == 0 || qi == 2) {
    const Real lon1 = 5*(M_PI/6), lat1 = 0, lon2 = -5*(M_PI/6), lat2 = 0;
    const Real zeta_c = 5e3, Z_t = 1e3, R_t = 0.5, lat_c = 0;
    
    const Real zeta = (zeta0 + zeta1)/2;
    Real lon, lat, d1, d2;
    xyz2lonlat(m.c.x[ci], m.c.y[ci], m.c.z[ci], lon, lat);
    const auto zdist2 = square((zeta - zeta_c)/Z_t);
    const auto dist_fn = [&] (const Real loni, const Real lati) {
      return std::min(1.0, (square(great_circle_dist(lon, lat, loni, lati)/R_t) +
                            zdist2));
    };
    d1 = dist_fn(lon1, lat1);
    d2 = dist_fn(lon2, lat2);
    assert(d1 >= 0 && d1 <= 1);
    assert(d2 >= 0 && d2 <= 1);

    if (qi == 0)
      q = (2 + std::cos(M_PI*d1) + std::cos(M_PI*d2))/2;
    else {
      assert(qi == 2);
      q = 0.1;
      if (d1 < 0.5 || d2 < 0.5) q = 1;
      if (zeta > zeta_c && std::abs(lat - lat_c) < 1.0/8) q = 0.1;            
    }
  } else {
    if (qi == 1)
      q = 0.9 - 0.8*square(get_ic_dcmip12t11(m, 0, ci, zeta0, zeta1));
    else if (qi == 4)
      q = 1;
    else {
      assert(qi == 3);
      q = 0;
      for (int i = 0; i < 3; ++i)
        q += get_ic_dcmip12t11(m, i, ci, zeta0, zeta1);
      q = 1 - 0.3*q;
    }
  }

  assert(q >= 0 && q <= 1);
  return q;
}

static Real get_max_cell_area (const MpasMesh& mm, Comm& he,
                               const Real* cell_area) {
  if ( ! cell_area) cell_area = mm.g.areaCell;
  Real lmax = 0, gmax;
  for (int i = 0; i < mm.ncown; ++i) lmax = std::max(lmax, cell_area[i]);
  he.all_reduce(Comm::MAX, 1, &lmax, &gmax);
  return gmax;
}

static Int get_glbl_num_cells (const MpasMesh& mm, Comm& he) {
  Int gncell;
  he.all_reduce(Comm::SUM, 1, &mm.ncown, &gncell);
  return gncell;
}

template <typename A3>
static Real get_max_mag (Comm& he, const A3& a, const int qi, const int nlvl,
                         const int ncell) {
  Real mv = 0;
  for (int ci = 0; ci < ncell; ++ci)
      for (int li = 0; li < nlvl; ++li)
        mv = std::max(mv, std::abs(a(qi,li,ci)));
  Real gmv;
  he.all_reduce(Comm::MAX, 1, &mv, &gmv);
  return gmv;
}

template <typename Accum> static void
get_masses (const MpasMesh& mm, Comm& he, const Int qi, const Int nlvl,
            const Int ncown, const Real* cell_area, const Ar3& Qt, const Ar3& Q,
            Accum& mass_beg, Accum& mass_end) {
  for (int ci = 0; ci < ncown; ++ci) {
    const Real* zint = mm.l.get_interface_z(ci);
    for (int k = 0; k < nlvl; ++k) {
      Real dz = zint ? zint[k+1] - zint[k] : 1;
      if (dz == 0) dz = 1; // 2D case.
      mass_beg += dz*cell_area[ci]*Qt(qi,k,ci);
      mass_end += dz*cell_area[ci]*Q (qi,k,ci);
    }
  }
}

void cmp (const MpasMesh& mm, Comm& he, const Int qi, const Int lvl, const Int nlvl,
          const Int ncown, const Ar2* rhot, const Ar3& Qt, const Ar2* rho, const Ar3& Q,
          CmpOutput& out, const Real* cell_area) {
  typedef order_invariant_accum::AccumType<double,true,false>::type Accum;
  assert( ! (rhot || rho) || (rhot && rho));
  if ( ! cell_area) cell_area = mm.g.areaCell;

  const auto ncell = get_glbl_num_cells(mm, he);
  const auto max_area = get_max_cell_area(mm, he, cell_area);
  Real max_scalar, maxmag; {
    const auto get_scalar = [&] (const int qi, const int li, const int ci) {
      const Real
        qt = rhot ? Qt(qi,li,ci)/(*rhot)(li,ci) : Qt(qi,li,ci),
        q = rho ? Q(qi,li,ci)/(*rho)(li,ci) : Q(qi,li,ci);
      return std::max(std::abs(qt), std::abs(q));
    };
    max_scalar = get_max_mag(he, get_scalar, qi, nlvl, ncown);
    maxmag = ncell*max_area*(lvl == -1 ? nlvl : 1)*max_scalar;
  }

  // Errors.
  Accum l2_num(maxmag), l2_den(maxmag), l1_num(maxmag), l1_den(maxmag);
  Real linf_num = -1e30, linf_den = -1e30;
  Real lo = 1e30, hi = -1e30, lot = lo, hit = hi;
  const auto f = [&] (const int ci, const int li) {
    const Real
      qt = rhot ? Qt(qi,li,ci)/(*rhot)(li,ci) : Qt(qi,li,ci),
      q = rho ? Q(qi,li,ci)/(*rho)(li,ci) : Q(qi,li,ci);
    l1_num += cell_area[ci]*std::abs(q - qt);
    l1_den += cell_area[ci]*std::abs(qt);
    l2_num += cell_area[ci]*square(q - qt);
    l2_den += cell_area[ci]*square(qt);
    linf_num = std::max(linf_num, std::abs(q - qt));
    linf_den = std::max(linf_den, std::abs(qt));
    if (rho) {
      lot = std::min(lot, qt);
      hit = std::max(hit, qt);
      lo = std::min(lo, q);
      hi = std::max(hi, q);
    }
  };
  if (lvl == -1) { // 3D
    for (int ci = 0; ci < ncown; ++ci)
      for (int li = 0; li < nlvl; ++li)
        f(ci, li);
  } else { // 2D
    for (int ci = 0; ci < ncown; ++ci)
      f(ci, lvl);
  }

  // Mass is always treated 3D b/c the CEDR acts on the 3D problem.
  Real mass_maxmag; {
    Real mass_beg = 0, mass_end = 0;
    get_masses(mm, he, qi, nlvl, ncown, cell_area, Qt, Q, mass_beg, mass_end);
    Real l[] = {mass_beg, mass_end}, g[2];
    he.all_reduce(Comm::SUM, 2, l, g);
    mass_maxmag = 2*std::max(g[0], g[1]);
  }
  Accum mass_beg(mass_maxmag), mass_end(mass_maxmag);
  get_masses(mm, he, qi, nlvl, ncown, cell_area, Qt, Q, mass_beg, mass_end);

  // Reductions.
  {
    Accum::IRep l[6], g[6];
    l[0] = order_invariant_accum::get_reduction_representation(mass_end);
    l[1] = order_invariant_accum::get_reduction_representation(mass_beg);
    l[2] = order_invariant_accum::get_reduction_representation(l1_num);
    l[3] = order_invariant_accum::get_reduction_representation(l1_den);
    l[4] = order_invariant_accum::get_reduction_representation(l2_num);
    l[5] = order_invariant_accum::get_reduction_representation(l2_den);
    he.all_reduce(Comm::SUM, 6, l, g);
    {
      const auto e = order_invariant_accum::make_accum_from_existing(mass_end, g[0]);
      const auto b = order_invariant_accum::make_accum_from_existing(mass_beg, g[1]);
      out.mass = (e() - b())/b();
    }
    {
      const auto n = order_invariant_accum::make_accum_from_existing(l1_num, g[2]);
      const auto d = order_invariant_accum::make_accum_from_existing(l1_den, g[3]);      
      out.l1 = n()/d();
    }
    {
      const auto n = order_invariant_accum::make_accum_from_existing(l2_num, g[4]);
      const auto d = order_invariant_accum::make_accum_from_existing(l2_den, g[5]);
      out.l2 = std::sqrt(n()/d());
    }
  }
  {
    Real vl[] = {lo, lot, hi, hit, linf_num, linf_den}, vg[6];
    he.all_reduce(Comm::MIN, 2, vl, vg);
    he.all_reduce(Comm::MAX, 4, vl+2, vg+2);
    out.lo = std::max(0.0, vg[1] - vg[0]);
    out.hi = std::max(0.0, vg[2] - vg[3]);
    out.linf = vg[4]/vg[5];
  }
}

Real VelocityField::max_vel () const {
  // As usual, assume radius 1.
  return 1.25*(2*M_PI)/cycle_time();
}

Real VelocityField::mlvl_to_z (const int nmlvl, const int mlvl) const {
  return 0.5*(ilvl_to_z(nmlvl+1, mlvl) + ilvl_to_z(nmlvl+1, mlvl+1));
}

Real VelocityField::get_ic_shape_fn (const Real z0, const Real z1) const {
  if (is_3d()) {
    const Real
      zetab = ilvl_to_z(2, 0),
      zetat = ilvl_to_z(2, 1),
      zeta  = (z0 + z1)/2,
      c     = (zeta - zetab)/(zetat - zetab);
//#define SHAPE_FN_STAIR_STEP
#ifndef SHAPE_FN_STAIR_STEP
    return 0.5*(1 - std::cos(2*M_PI*c));
#else
    // Make vertical structure clear by using a stair-step shape fn.
    Real f = 0.5*(1 - std::cos(2*M_PI*c));
    return std::round(1 + 4*f)/5;
#endif
  } else {
    return 1;
  }
}

Rotation::Rotation (const bool zaxis) {
  Omega[0] = 0.2; Omega[1] = 0.7; Omega[2] = 1;
  if (zaxis) { Omega[0] = Omega[1] = 0; }
  geo::normalize(Omega);
  const Real a = (2*M_PI)/(12*86400);
  for (int i = 0; i < 3; ++i) Omega[i] *= a;
}

void Rotation
::operator() (const Real t, const Real lon, const Real lat, Real, Int,
              Real& u, Real& v, Real& w) const {
  Real x, y, z;
  lonlat2xyz(lon, lat, x, y, z);
  const Vec r = {x,y,z};
  Vec vel; geo::cross(Omega, r, vel);
  const Vec axis = {0,0,1};
  Vec east; geo::cross(axis, r, east); geo::normalize(east);
  Vec north; geo::cross(r, (const Real*) east, north); geo::normalize(north);
  u = geo::dot(vel, east);
  v = geo::dot(vel, north);
  w = 0;
}

void NonDivergentWindField
::operator() (const Real t, const Real lon, const Real lat, Real, Int,
              Real& u, Real& v, Real& w) const {
  const Real theta = lat, lambda = lon;
  const Real
    T = 12*86400,
    R = 1,
    lambda_p = lambda - 2*M_PI*t/T,
    costh = std::cos(theta),
    cost = std::cos(M_PI*t/T);
  u = R/T*(10*square(std::sin(lambda_p))*std::sin(2*theta)*cost + 2*M_PI*costh);
  v = 10*R/T*std::sin(2*lambda_p)*costh*cost;
  w = 0;
}

void DivergentWindField
::operator() (const Real t, const Real lon, const Real lat, Real, Int,
              Real& u, Real& v, Real& w) const {
  const Real theta = lat, lambda = lon;
  const Real
    T = 12*86400,
    R = 1,
    lambda_p = lambda - 2*M_PI*t/T,
    costh = std::cos(theta),
    cost = std::cos(M_PI*t/T);
  u = R/T*(-5*square(std::sin(0.5*lambda_p))*std::sin(2*theta)*square(costh)*cost
           + 2*M_PI*costh);
  v = 2.5*R/T*std::sin(lambda_p)*cube(costh)*cost;
  w = 0;
}

void Multiple2DWindFields
::operator() (const Real t, const Real lon, const Real lat, Real, const Int level,
              Real& u, Real& v, Real& w) const {
  switch (std::abs(level) % 3) {
    case 0: { NonDivergentWindField f; f(t, lon, lat, 0, 0, u, v, w); }; break;
    case 1: {    DivergentWindField f; f(t, lon, lat, 0, 0, u, v, w); }; break;
    case 2: {              Rotation f; f(t, lon, lat, 0, 0, u, v, w); }; break;
  }
}

struct Data3D {
  static constexpr Real
    ztop = 12000,
    T0 = 300,
    Rd = 287.04,
    g  = 9.80616,
    p0 = 100000,
    H  = Rd*T0/g;

  static Real z2p (const Real z) { return p0*std::exp(-z/H); }
  static Real get_rho (const Real p) { return p/(Rd*T0); }
  static Real get_ptop () { return z2p(ztop); }
};

static void
calc_nondiv (const Real bs, const Real strange_two, const Real my_two,
             const Real t, const Real lon, const Real lat, Real z,
             Real& u, Real& v, Real& w) {
  typedef Data3D d;
  
#ifndef NDEBUG
  static const Real tol = 1e4*std::numeric_limits<Real>::epsilon();
  assert(z >= -tol && z <= (1+tol)*d::ztop);
#endif
  if (z < 0) z = 0;
  if (z > d::ztop) z = d::ztop;

  const Real
    R  = 1,
    T  = 12*86400,
    u0 = 2*M_PI*R/T,
    k0 = 10*R/T,
    // I put in a 2 to make the trajectories more exciting.
    omega0 = my_two*(23000*M_PI)/T,
    theta = lat, lambda = lon,
    lambda_p = lambda - 2*M_PI*t/T,
    p       = d::z2p(z),
    ptop    = d::get_ptop(),
    costh   = std::cos(theta),
    cost    = std::cos(M_PI*t/T),
    // This has a 2 in it in the doc, but integrating the trajectories in time
    // shows this then doesn't return a point to its original position at day
    // 12, as is intended. It's also clear that with the 2, the ud and w fields
    // in isolation don't reverse. Run `hy mpasdbg.hy vis` to see trajectories.
    cos2t   = std::cos(strange_two*M_PI*t/T),
    sinlamp = std::sin(lambda_p);

  const Real
    s = (1.0 + std::exp((ptop-d::p0)/(bs*ptop)) -
         std::exp((p-d::p0)/(bs*ptop)) - std::exp((ptop-p)/(bs*ptop))),
    s_p = (-std::exp((p-d::p0)/(bs*ptop)) +
           std::exp((ptop-p)/(bs*ptop)))/(bs*ptop),
    ud = (omega0*R)*std::cos(lambda_p)*square(costh)*cos2t*s_p;

#if 1
  u = k0*square(sinlamp)*std::sin(2*theta)*cost + u0*costh + ud;
  v = k0*std::sin(2*lambda_p)*costh*cost;
  w = (-1/(d::g*d::get_rho(p)))*omega0*sinlamp*costh*cos2t*s;
#else
  // Use this to convergence-test the vertical direction in isolation.
  u = v = 0;
  w = (-1/(d::g*d::get_rho(p)))*omega0*cos2t*s;
#endif  
}

Real Dcmip20212Test11WindField
::ilvl_to_z (const int nilvl, const int ilvl) const {
  return (Real(ilvl)/std::max(1,nilvl-1))*Data3D::ztop;
}

Real Dcmip20212Test11WindField::get_rho (const Real z0, const Real z1) const {
  return Data3D::get_rho(Data3D::z2p((z0 + z1)/2));
}

Real Dcmip20212Test11WindField
::get_ic_shape_fn (const Real z0, const Real z1) const { return 1; }

void Dcmip20212Test11WindField
::operator() (const Real t, const Real lon, const Real lat, Real z,
              Int, Real& u, Real& v, Real& w) const {
  calc_nondiv(0.2, 2, 1, t, lon, lat, z, u, v, w);
}

Real NonDivergent3DWindField::ilvl_to_z (const int nilvl, const int ilvl) const {
  return (Real(ilvl)/std::max(1,nilvl-1))*Data3D::ztop;
}

Real NonDivergent3DWindField::get_rho (const Real z0, const Real z1) const {
  return Data3D::get_rho(Data3D::z2p((z0 + z1)/2));
}

void NonDivergent3DWindField
::operator() (const Real t, const Real lon, const Real lat, Real z,
              Int, Real& u, Real& v, Real& w) const {
  calc_nondiv(1, 1, 2, t, lon, lat, z, u, v, w);
}

struct Data3DVertP {
  typedef Data3D d;
  static constexpr auto p0 = d::p0;
  static constexpr auto ztop = d::ztop;
  // Vertical coord to p. Since the v coord is already p, just return it.
  static Real z2p (const Real z) { return z; }
  static Real p2height (const Real p) { return -d::H*std::log(p/d::p0); }
  // This is pseudodensity, dz/dp rho \propto 1.
  static Real get_rho (const Real) { return 1; }
  static Real get_ptop () { return d::get_ptop(); }
};

static void
calc_nondiv_vertp (const Real bs, const Real strange_two, const Real my_two,
                   const Real t, const Real lon, const Real lat, Real z,
                   Real& u, Real& v, Real& w) {
  typedef Data3DVertP d;

  const Real ptop = d::get_ptop();
  
#ifndef NDEBUG
  static const Real tol = 1e4*std::numeric_limits<Real>::epsilon();
  assert(z >= (1 - tol)*ptop && z <= (1+tol)*d::p0);
#endif
  if (z < ptop) z = ptop;
  if (z > d::p0) z = d::p0;

  const Real
    R  = 1,
    T  = 12*86400,
    u0 = 2*M_PI*R/T,
    k0 = 10*R/T,
    omega0 = my_two*(23000*M_PI)/T,
    theta = lat, lambda = lon,
    lambda_p = lambda - 2*M_PI*t/T,
    p       = d::z2p(z),
    costh   = std::cos(theta),
    cost    = std::cos(M_PI*t/T),
    cos2t   = std::cos(strange_two*M_PI*t/T),
    sinlamp = std::sin(lambda_p);

  const Real
    s = (1.0 + std::exp((ptop-d::p0)/(bs*ptop)) -
         std::exp((p-d::p0)/(bs*ptop)) - std::exp((ptop-p)/(bs*ptop))),
    s_p = (-std::exp((p-d::p0)/(bs*ptop)) +
           std::exp((ptop-p)/(bs*ptop)))/(bs*ptop),
    ud = (omega0*R)*std::cos(lambda_p)*square(costh)*cos2t*s_p;

  u = k0*square(sinlamp)*std::sin(2*theta)*cost + u0*costh + ud;
  v = k0*std::sin(2*lambda_p)*costh*cost;
  w = omega0*sinlamp*costh*cos2t*s;
}

NonDivergent3DVertPWindField::NonDivergent3DVertPWindField () {
#if 0
  const int nilvl = 10;
  for (int i = 0; i < nilvl; ++i) {
    const auto ivert = ilvl_to_z(nilvl, i);
    const auto height = Data3DVertP::p2height(ivert);
    Real dp = 0, dheight = 0;
    if (i < nilvl-1) {
      const auto ip1vert = ilvl_to_z(nilvl, i+1);
      dp = ip1vert - ivert;
      dheight = height - Data3DVertP::p2height(ip1vert);
    }
    pr(puf(i) pu(ivert) pu(height) pu(dp) pu(dheight));
  }
#endif
}

Real NonDivergent3DVertPWindField::ilvl_to_z (const int nilvl, const int ilvl) const {
  const auto z = (Real(nilvl-1-ilvl)/std::max(1,nilvl-1))*Data3DVertP::ztop;
  // Data3D, not Data3DVertP, is indeed intended.
  return Data3D::z2p(z);
}

Real NonDivergent3DVertPWindField::get_rho (const Real z0, const Real z1) const {
  return Data3DVertP::get_rho(Data3DVertP::z2p((z0 + z1)/2));
}

Real NonDivergent3DVertPWindField
::get_ic_shape_fn (const Real p0, const Real p1) const {
  const Real
    zetab = 0,
    zetat = Data3DVertP::ztop,
    z0    = Data3DVertP::p2height(p1),
    z1    = Data3DVertP::p2height(p0),
    zeta  = (z0 + z1)/2,
    c     = (zeta - zetab)/(zetat - zetab);
#ifndef SHAPE_FN_STAIR_STEP
  return 0.5*(1 - std::cos(2*M_PI*c));
#else
  Real f = 0.5*(1 - std::cos(2*M_PI*c));
  return std::round(1 + 4*f)/5;
#endif
}

void NonDivergent3DVertPWindField
::operator() (const Real t, const Real lon, const Real lat, Real z,
              Int, Real& u, Real& v, Real& w) const {
  calc_nondiv_vertp(1, 1, 2, t, lon, lat, z, u, v, w);
}

Real TestWindField::ilvl_to_z (const int nlvl, const int lvl) const {
  assert(lvl >= 0 && lvl < nlvl);
  return ((lvl+0.5)/nlvl)*Data3D::ztop;
}

// Study the 3D field in a layer near the bottom boundary.
void TestWindField
::operator() (const Real t, const Real lon, const Real lat, Real,
              Int, Real& u, Real& v, Real& w) const {
  typedef Data3D d;
  const Real z = 0.05*d::ztop;
    
  static constexpr Real
    R  = 1,
    T  = 12*86400,
    u0 = 2*M_PI*R/T,
    k0 = 10*R/T,
    omega0 = (2*23000*M_PI)/T,
    bs = 1;

  const Real
    theta = lat, lambda = lon,
    lambda_p = lambda - 2*M_PI*t/T,
    p    = d::z2p(z),
    ptop = d::get_ptop(),
    costh   = std::cos(theta),
    cost    = std::cos(M_PI*t/T),
    sinlamp = std::sin(lambda_p);

  const Real
    s_p = (-std::exp((p-d::p0)/(bs*ptop)) +
           std::exp((ptop-p)/(bs*ptop)))/(bs*ptop),
    ud = (omega0*R)*std::cos(lambda_p)*square(costh)*cost*s_p;

  u = k0*square(sinlamp)*std::sin(2*theta)*cost + u0*costh + ud;
  v = k0*std::sin(2*lambda_p)*costh*cost;
  w = 0;
}

Real HadleyWindField::ilvl_to_z (const int nilvl, const int ilvl) const {
  return (Real(ilvl)/std::max(1,nilvl-1))*Data3D::ztop;
}

Real HadleyWindField::get_ic_shape_fn (const Real zlo, const Real zhi) const {
  const Real z1 = 2e3, z2 = 5e3, z0 = 0.5*(z1 + z2);
#if 1
  // Use the zeta midpoint value.
  const Real z = (zlo + zhi)/2;
  return ((z > z1 && z < z2) ?
          0.5*(1 + std::cos(2*M_PI*((z - z0)/(z2 - z1)))) :
          0);
#elif 0
  assert(0);
  // Use the average value:
  //   int[ (1 + cos(2 pi (z - z0)/(z2 - z1)))/2; zlo, zhi ]/(zhi - zlo)
  //   diffat[ z/2 + (z2 - z1)/(4 pi) sin(2 pi (z - z0)/(z2 - z1));
  //           zlo, zhi ]/(zhi - zlo)
  if (zlo >= z2 || zhi <= z1) return 0;
  const auto f = [&] (const Real z) {
    return 0.5*(z + (z2 - z1)/(2*M_PI)*std::sin(2*M_PI*((z - z0)/(z2 - z1))));
  };
  const Real zlo1 = std::max(zlo, z1), zhi1 = std::min(zhi, z2);
  return (f(zhi1) - f(zlo1))/(zhi - zlo);
#else
  assert(0);
  // Gaussian.
  const Real z = (zlo + zhi)/2;
  // exp(-b (z1-z0)^2) = 0.01
  // => b = -log(1/2)/(z1-z0)^2
  const Real b = -log(0.1)/square(z1-z0);
  return exp(-b*square(z - z0));
#endif
}

Real HadleyWindField::get_rho (const Real z0, const Real z1) const {
  return Data3D::get_rho(Data3D::z2p((z0 + z1)/2));
}

void HadleyWindField
::operator() (const Real t, const Real lon, const Real lat, Real z,
              Int, Real& u, Real& v, Real& w) const {
  typedef Data3D d;

  const Real
    T = cycle_time(),
    aref = 6.37122e6,
    u0 = 40,
    w0 = 0.15,
    K = 5,
    theta = lat,
    rho0 = d::get_rho(d::p0),
    rho = d::get_rho(d::z2p(z)),
    costh = std::cos(theta), sinth = std::sin(theta),
    sinKth = std::sin(K*theta), cosKth = std::cos(K*theta),
    cost = std::cos(M_PI*t/T);

  u = (u0/aref)*costh;
  v = (-((w0*M_PI*rho0)/(K*d::ztop*rho))*
       costh*sinKth*std::cos(M_PI*z/d::ztop)*cost);
  w = (((w0*rho0)/(K*rho))*
       (-2*sinKth*sinth + K*costh*cosKth)*std::sin(M_PI*z/d::ztop)*cost);
}

VelocityField::Ptr vf_from_short (char sym, const bool verbose) {
  switch (sym) {
  case 'n': return std::make_shared<NonDivergentWindField>();
  case 'd': return std::make_shared<DivergentWindField>();
  case 'r': return std::make_shared<Rotation>();
  case 'a': return std::make_shared<Multiple2DWindFields>();
  case '3': return std::make_shared<NonDivergent3DWindField>();
  case '4': return std::make_shared<NonDivergent3DVertPWindField>();
  case '1': return std::make_shared<Dcmip20212Test11WindField>();
  case 't': return std::make_shared<TestWindField>();
  case 'h': return std::make_shared<HadleyWindField>();
  default:
    if (verbose) printf("flow field not supported: %c\n", sym);
    assert(0);
    return nullptr;
  }
}

InitialCondition::Shape ic_from_short (char sym, const bool verbose) {
  switch (sym) {
  case 'g': return InitialCondition::Shape::GaussianHills;
  case 'c': return InitialCondition::Shape::CosineBells;
  case 's': return InitialCondition::Shape::SlottedCylinders;
  case 't': return InitialCondition::Shape::XYZTrig;
  case '1': return InitialCondition::Shape::Constant;
  default:
    if (verbose) printf("slot 1: ic not supported: %c\n", sym);
    assert(0);
    return InitialCondition::Shape::XYZTrig;
  }
}

Int calc_default_nstep (const MpasMesh& m, Comm& he, const VelocityField& vf) {
  Real R, D;
  get_avg_radius_and_cell_spacing(m, he, R, D);
  const Real T = vf.cycle_time(), max_vel = vf.max_vel(), cfl = (D/R)/max_vel;
  static const Int hop = 5; // ~nbr cells to go in a time step
  return std::max<Real>(
    10,
    ((T/cfl)/
     (he.parallel() || he.szhalo() > 0 ?
      // If parallel, reduce hop if szhalo isn't large enough. Even if serial,
      // if szhalo is > 0, use it so that serial and parallel runs get the same
      // settings.
      std::min(std::max(1, he.szhalo()-2), hop):
      hop) + 1));
}

// One of nhalo or nsubstep should be < 0 but not both.
static
void calc_substep_datum (const MpasMesh& m, Comm& he, const VelocityField& vf,
                         const int nstep, int& nhalo, int& nsubstep) {
  assert(nhalo < 0 || nsubstep < 0);
  assert(nhalo > 0 || nsubstep > 0);
  Real R, D;
  get_avg_radius_and_cell_spacing(m, he, R, D);
  const Real
    max_vel = vf.max_vel(),
    cfl = (D/R)/max_vel;
  const int n_ff_step = (int) (vf.cycle_time()/(0.4*cfl) + 1);
  const int n = nhalo < 0 ? nsubstep : nhalo;
  const int npart = (n_ff_step + nstep*n - 1)/(nstep*n);
  if (nhalo < 0) nhalo = npart;
  else nsubstep = npart;
}

void calc_substep_data (const MpasMesh& m, Comm& he, const VelocityField& vf,
                        const int nstep, int& traj_nhalo, int& nsubstep) {
  if (nsubstep < 0) {
    if (traj_nhalo < 1) traj_nhalo = max_nhalo;
    nsubstep = -1;
    calc_substep_datum(m, he, vf, nstep, traj_nhalo, nsubstep);
  } else {
    int nhalo = -1;
    calc_substep_datum(m, he, vf, nstep, nhalo, nsubstep);
    if (nhalo > max_nhalo) {
      traj_nhalo = max_nhalo;
      nsubstep = -1;
      calc_substep_datum(m, he, vf, nstep, traj_nhalo, nsubstep);
    } else
      traj_nhalo = std::max(traj_nhalo, nhalo);
  }
}

void set_interface_z_all_cols (const VelocityField& vf, MpasMesh& mm) {
  auto& l = mm.l;
  std::vector<Real> zint(l.nmax+1);
  for (int i = 0; i <= l.nmax; ++i) zint[i] = vf.ilvl_to_z(l.nmax+1, i);
  l.set_interface_z_all_cols(zint.data());
}

Real get_level_shape_fn (const Real x, const Real ctr, const Real a) {
  assert(x >= -1 && x <= 1);
  assert(ctr >= -1 && ctr <= 1);
  assert(a > 0);
  const auto g = [&] (const Real x) { return 1/(1 + std::exp(-a*x)); };
  const Real c = -ctr, gcm1 = g(c-1);
  return (g(c+x) - gcm1)/(g(c+1) - gcm1);
}
