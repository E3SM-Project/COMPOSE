// mycpp -g test.cpp; if [ $? == 0 ]; then ./a.out -m | grep "mat=1" > foo.m; fi
// >> clf;msik('draw_test_output','foo');ic

// ko=/home/ambradl/lib/kokkos/cpu; mycpp -DSIQK_USE_KOKKOS -I$ko/include -L$ko/lib -fopenmp test.cpp -lkokkos -ldl

#ifdef SIQK_USE_KOKKOS
# include "Array_Kokkos.hpp"
#else
# include "Array_raw.hpp"
#endif
#include "sik.hpp"
using namespace siqk;

typedef SphereGeometry Geo;

template <typename T>
void copy (Array2D<T>& d, const Array2D<const T>& s) {
  for (auto i = zero(s); i < s.n(); ++i)
    for (auto j = zero(s); j < s.m(); ++j)
      d(j,i) = s(j,i);
}

static void
write_matlab (const std::string& name, const Array2D<const double>& p) {
  printf("mat=1; %s = [", name.c_str());
  for (int ip = zero(p); ip < p.n(); ++ip)
    printf(" %1.15e %1.15e %1.15e;", p(0,ip), p(1,ip), p(2,ip));
  printf("].';\n");
}

static void
write_matlab (const std::string& name, const Array2D<const double>& p,
              const Array2D<const int>& e) {
  printf("mat=1; %s.p = [", name.c_str());
  for (int ip = zero(p); ip < p.n(); ++ip)
    printf(" %1.15e %1.15e %1.15e;", p(0,ip), p(1,ip), p(2,ip));
  printf("].';\n");
  printf("mat=1; %s.n = [", name.c_str());
  for (int ie = zero(e); ie < e.n(); ++ie)
    printf(" %d %d %d %d;", e(0,ie)+1, e(1,ie)+1, e(2,ie)+1, e(3,ie)+1);
  printf("].';\n");
}

static void make_planar_mesh (Array2D<double>& p, Array2D<int>& e,
                              const int n) {
  const double d = std::sqrt(0.5);
  e.reset(4, n*n);
  p.reset(3, (n+1)*(n+1));
  p.set(0);
  for (int iy = 0; iy < n+1; ++iy)
    for (int ix = 0; ix < n+1; ++ix) {
      const auto idx = (n+1)*iy + ix;
      p(0,idx) = 2*(static_cast<double>(ix)/n - 0.5)*d;
      p(1,idx) = 2*(static_cast<double>(iy)/n - 0.5)*d;
    }
  for (int iy = 0; iy < n; ++iy)
    for (int ix = 0; ix < n; ++ix) {
      const auto idx = n*iy + ix;
      e(0,idx) = (n+1)*iy + ix;
      e(1,idx) = (n+1)*iy + ix+1;
      e(2,idx) = (n+1)*(iy+1) + ix+1;
      e(3,idx) = (n+1)*(iy+1) + ix;
    }
}

static void project_onto_sphere (Array2D<double>& p) {
  for (auto ip = zero(p); ip < p.n(); ++ip) {
    p(2,ip) = 1;
    SphereGeometry::normalize(p(ip));
  }
}

static void
perturb_mesh (Array2D<double>& p, Array2D<int>& e, const double angle,
              const double xlate, const double ylate) {
  const double cr = std::cos(angle), sr = std::sin(angle);
  for (auto ip = zero(p); ip < p.n(); ++ip) {
    const double x = p(0,ip), y = p(1,ip);
    p(0,ip) =  cr*x - sr*y + xlate;
    p(1,ip) = -sr*x + cr*y + ylate;
  }  
}

static void fill_quad (const Array2D<const double>& p, Array2D<double>& poly) {
  const int n = static_cast<int>(std::sqrt(p.n() - 1));
  copy(poly(0), p(0), 3);
  copy(poly(1), p(n), 3);
  copy(poly(2), p(p.n() - 1), 3);
  copy(poly(3), p(p.n() - 1 - n), 3);
}

// Area of the outline of (p,e) clipped against the outline of (cp,ce).
static double
calc_true_area (const Array2D<const double>& cp, const Array2D<const int>& ce,
                const Array2D<const double>& p, const Array2D<const int>& e,
                const bool wm) {
  Array2D<double> clip_poly(3, 4), poly(3, 4), nml(3, 4);
  fill_quad(cp, clip_poly);
  fill_quad(p, poly);
  for (int i = 0; i < 4; ++i)
    Geo::edge_normal(clip_poly(i), clip_poly((i+1) % 4), nml(i));
  Array2D<double> vo(3, sh::max_nvert);
  int no;
  sh::clip_against_poly<Geo>(clip_poly, nml, poly, 4, vo, no);
  Array2D<const double> intersection(3, no, vo.data());
  if (wm) {
    write_matlab("clip_poly", clip_poly);
    write_matlab("poly", poly);
    write_matlab("intersection", intersection);
  }
  return Geo::calc_area(intersection);
}

static int
run (const int n, const double angle, const double xlate, const double ylate,
     const bool wm) {
  // Make the clip mesh.
  Array2D<double> cp;
  Array2D<int> ce;
  make_planar_mesh(cp, ce, n);

  // Make a perturbed mesh.
  Array2D<double> p(cp.m(), cp.n());
  Array2D<int> e(ce.m(), ce.n());
  copy<double>(p, cp);
  copy<int>(e, ce);
  perturb_mesh(p, e, angle, xlate, ylate);

  // Project these meshes onto the sphere.
  project_onto_sphere(cp);
  project_onto_sphere(p);

  // True intersection area from quadrilateral boundary of the mesh.
  const double ta = calc_true_area(cp, ce, p, e, wm);
  // Area from the sum over the common refinement polygons.
  const double a = test::test_area_ot<Geo>(cp, ce, p, e);

  // Report information.
  const double re = std::abs(a - ta)/ta;
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
  int n;
  double angle, xlate, ylate;
  bool write_matlab;

  Input (int argc, char** argv)
    : n(5), angle(M_PI*1e-1), xlate(1e-1), ylate(1e-1), write_matlab(false)
  {
    for (int i = 1; i < argc; ++i) {
      const std::string& token = argv[i];
      if (eq(token, "-n"))
        n = atoi(argv[++i]);
      if (eq(token, "-m", "--write-matlab"))
        write_matlab = true;
    }

    print(std::cout);
  }

  void print (std::ostream& os) {
    os << "n (-n): " << n
       << "\n";
  }
};

int main (int argc, char** argv) {
#ifdef SIQK_USE_KOKKOS
  Kokkos::initialize(argc, argv);
#endif
  {
    Input in(argc, argv);
    run(in.n, in.angle, in.xlate, in.ylate, in.write_matlab);  
  }
#ifdef SIQK_USE_KOKKOS
  Kokkos::finalize_all();
#endif
}
