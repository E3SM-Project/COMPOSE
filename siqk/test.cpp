// ko=/home/ambradl/lib/kokkos/cpu; mycpp -I$ko/include -L$ko/lib -fopenmp test.cpp -lkokkos -ldl -Wall -pedantic
// ./a.out -m | grep "mat=1" > foo.m
// >> msik('draw_test_output', 'foo');

#include "siqk_intersect.hpp"
using namespace siqk;

template <typename CV3s>
static void
write_matlab (const std::string& name, const CV3s& p) {
  printf("mat=1; %s = [", name.c_str());
  for (Int ip = 0; ip < nslices(p); ++ip)
    printf(" %1.15e %1.15e %1.15e;", p(ip,0), p(ip,1), p(ip,2));
  printf("].';\n");
}

template <typename CV3s, typename CIs>
static void
write_matlab (const std::string& name, const CV3s& p, const CIs& e) {
  printf("mat=1; %s.p = [", name.c_str());
  for (Int ip = 0; ip < nslices(p); ++ip)
    printf(" %1.15e %1.15e %1.15e;", p(ip,0), p(ip,1), p(ip,2));
  printf("].';\n");
  printf("mat=1; %s.n = [", name.c_str());
  for (Int ie = 0; ie < nslices(e); ++ie)
    printf(" %d %d %d %d;", e(ie,0)+1, e(ie,1)+1, e(ie,2)+1, e(ie,3)+1);
  printf("].';\n");
}

static void make_planar_mesh (Vec3s::HostMirror& p, Idxs::HostMirror& e,
                              const Int n) {
  const Real d = std::sqrt(0.5);
  ko::resize(e, n*n, 4);
  ko::resize(p, (n+1)*(n+1), 3);
  for (Int iy = 0; iy < n+1; ++iy)
    for (Int ix = 0; ix < n+1; ++ix) {
      const auto idx = (n+1)*iy + ix;
      p(idx,0) = 2*(static_cast<Real>(ix)/n - 0.5)*d;
      p(idx,1) = 2*(static_cast<Real>(iy)/n - 0.5)*d;
      p(idx,2) = 0;
    }
  for (Int iy = 0; iy < n; ++iy)
    for (Int ix = 0; ix < n; ++ix) {
      const auto idx = n*iy + ix;
      e(idx,0) = (n+1)*iy + ix;
      e(idx,1) = (n+1)*iy + ix+1;
      e(idx,2) = (n+1)*(iy+1) + ix+1;
      e(idx,3) = (n+1)*(iy+1) + ix;
    }
}

static void project_onto_sphere (Vec3s::HostMirror& p) {
  for (Int ip = 0; ip < nslices(p); ++ip) {
    p(ip,2) = 1;
    SphereGeometry::normalize(slice(p, ip));
  }
}

static void
perturb_mesh (Vec3s::HostMirror& p, Idxs::HostMirror& e,
              const Real angle, const Real xlate, const Real ylate) {
  const Real cr = std::cos(angle), sr = std::sin(angle);
  for (Int ip = 0; ip < nslices(p); ++ip) {
    const Real x = p(ip,0), y = p(ip,1);
    p(ip,0) =  cr*x - sr*y + xlate;
    p(ip,1) = -sr*x + cr*y + ylate;
  }  
}

static void fill_quad (const ConstVec3s::HostMirror& p,
                       Vec3s::HostMirror& poly) {
  const Int n = static_cast<int>(std::sqrt(nslices(p) - 1));
  copy(slice(poly, 0), slice(p, 0), 3);
  copy(slice(poly, 1), slice(p, n), 3);
  copy(slice(poly, 2), slice(p, nslices(p) - 1), 3);
  copy(slice(poly, 3), slice(p, nslices(p) - 1 - n), 3);
}

// Area of the outline of (p,e) clipped against the outline of (cp,ce).
template <typename Geo>
static Real calc_true_area (
  const ConstVec3s::HostMirror& cp, const ConstIdxs::HostMirror& ce,
  const ConstVec3s::HostMirror& p, const ConstIdxs::HostMirror& e,
  const bool wm)
{
  Vec3s::HostMirror clip_poly("clip_poly", 4, 3), poly("poly", 4, 3),
    nml("nml", 4, 3);
  fill_quad(cp, clip_poly);
  fill_quad(p, poly);
  for (Int i = 0; i < 4; ++i)
    Geo::edge_normal(slice(clip_poly, i), slice(clip_poly, (i+1) % 4),
                     slice(nml, i));
  Vec3s::HostMirror vo("vo", test::max_nvert, 3);
  Int no;
  {
    Vec3s::HostMirror wrk("wrk", test::max_nvert, 3);
    sh::clip_against_poly<Geo>(clip_poly, nml, poly, 4, vo, no, wrk);
  }
  if (wm) {
    write_matlab("clip_poly", clip_poly);
    write_matlab("poly", poly);
    write_matlab("intersection",
                 ko::subview(vo, std::pair<Int,Int>(0, no), ko::ALL()));
  }
  return Geo::calc_area_formula(vo, no);
}

template <typename Geo> void finalize_mesh (Vec3s::HostMirror& p) {}
template <> void finalize_mesh<SphereGeometry> (Vec3s::HostMirror& p) {
  project_onto_sphere(p);
}

template <typename Geo>
static int
run (const Int n, const Real angle, const Real xlate, const Real ylate,
     const bool wm) {
  Vec3s::HostMirror cp;
  Idxs::HostMirror ce;
  make_planar_mesh(cp, ce, n);

  Vec3s::HostMirror p("p", nslices(cp), szslice(cp));
  Idxs::HostMirror e("e", nslices(ce), szslice(ce));
  ko::deep_copy(p, cp);
  ko::deep_copy(e, ce);
  perturb_mesh(p, e, angle, xlate, ylate);

  finalize_mesh<Geo>(cp);
  finalize_mesh<Geo>(p);

  const Real ta = calc_true_area<Geo>(cp, ce, p, e, wm);
  const Real a = test::test_area_ot<Geo>(cp, ce, p, e);

  const Real re = std::abs(a - ta)/ta;
  fprintf(stderr, "true area %1.4e mesh area %1.4e relerr %1.4e\n", ta, a, re);
  if (wm) {
    write_matlab("cm", cp, ce);
    write_matlab("m", p, e);
  }
  return re < 1e-10 ? 0 : 1;
}

inline bool
eq (const std::string& a, const char* const b1, const char* const b2 = 0) {
  return (a == std::string(b1) || (b2 && a == std::string(b2)) ||
          a == std::string("-") + std::string(b1));
}

struct Input {
  Int n;
  Real angle, xlate, ylate;
  bool write_matlab, geo_sphere;

  Input (Int argc, char** argv)
    : n(5), angle(M_PI*1e-1), xlate(1e-1), ylate(1e-1), write_matlab(false),
      geo_sphere(true)
  {
    for (Int i = 1; i < argc; ++i) {
      const std::string& token = argv[i];
      if (eq(token, "-n")) n = atoi(argv[++i]);
      if (eq(token, "-m", "--write-matlab")) write_matlab = true;
      if (eq(token, "--plane")) geo_sphere = false;
      if (eq(token, "--xlate")) xlate = atof(argv[++i]);
      if (eq(token, "--ylate")) ylate = atof(argv[++i]);
      if (eq(token, "--angle")) angle = atof(argv[++i]);
    }

    print(std::cout);
  }

  void print (std::ostream& os) {
    os << "n (-n): " << n << "\n"
       << "write matlab (-m): " << write_matlab << "\n"
       << "planar geometry (--plane): " << ! geo_sphere << "\n"
       << "angle (--angle): " << angle << "\n"
       << "xlate (--xlate): " << xlate << "\n"
       << "ylate (--ylate): " << ylate << "\n";
  }
};

int main (int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  {
    Input in(argc, argv);
    Int nerr = 0;
    nerr += (in.geo_sphere ?
             run<SphereGeometry>(in.n, in.angle, in.xlate, in.ylate, in.write_matlab) :
             run<PlaneGeometry>(in.n, in.angle, in.xlate, in.ylate, in.write_matlab));
    std::cerr << (nerr ? "FAIL" : "PASS") << "ED\n";
  }
  Kokkos::finalize_all();
}
