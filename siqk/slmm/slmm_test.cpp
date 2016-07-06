#include "slmm_defs.hpp"
#include "slmm_mesh.hpp"
#include "slmm_gll.hpp"
#include "slmm_io.hpp"
#include "slmm_time_int.hpp"
#include "slmm_gallery.hpp"
#include "slmm_debug.hpp"
using namespace slmm;

struct Command {
  enum Enum {
    test_make_cubedsphere, test_gll, test_make_gll_mesh, test_time_int
  };
  static Enum from_string (const std::string& s) {
    if (s == "test_make_cubedsphere") return test_make_cubedsphere;
    if (s == "test_gll") return test_gll;
    if (s == "test_make_gll_mesh") return test_make_gll_mesh;
    if (s == "test_time_int") return test_time_int;
    throw std::runtime_error(s + " is not a command.");
  }
};

struct Input {
  Command::Enum command;
  Int n;
  Real angle;
  bool write_matlab, quiet;
  std::string fn_pre_out;

  Input(Int argc, char** argv);
  void print(std::ostream& os) const;
};

static Int test_make_cubedsphere (const Input& in) {
  const Int np = 4;
  Vec3s::HostMirror cp;
  Idxs::HostMirror c2n;
  mesh::make_cubedsphere(cp, c2n, in.n);
  Int nerr = 0;
  {
    const Int ne = mesh::impl::check_elem_normal_against_sphere(cp, c2n);
    if (ne) std::cerr << "FAIL: check_elem_normal_against_sphere\n";
    nerr += ne;
  }
  {
    Idxs::HostMirror c2e, e2n;
    mesh::impl::make_c2e_from_c2n(np, c2n, c2e, e2n);
    Int ne = 0;
    // Every edge has two cells, and each cell is a quad.
    if (nslices(e2n) != 4/2*nslices(c2n)) {
      ++ne;
      std::cerr << "FAIL: make_c2e_from_c2n\n";
    }
    nerr += ne;
  }
  if (in.write_matlab)
    write_matlab("cm", cp, c2n);
  return nerr;
}

static Int test_gll (const Input& in) {
  Int nerr = 0;
  const Real tol = 1e2*std::numeric_limits<Real>::epsilon();
  GLL gll;
  const Real* x, * wt;
  const Int np = 4;
  gll.get_coef(np, x, wt);
  Real sum = 0;
  for (Int i = 0; i < np; ++i)
    sum += wt[i];
  if (std::abs(2 - sum) > tol) ++nerr;
  for (Int j = 0; j < np; ++j) {
    Real gj[GLL::max_np]; gll.eval(np, x[j], gj);
    for (Int i = 0; i < np; ++i) {
      if (j == i) continue;
      if (std::abs(gj[i]) > tol) ++nerr;
    }
  }
  return nerr;
}

static Int test_make_gll_mesh (const Input& in) {
  const Int np = 4;
  Vec3s::HostMirror geo_p, cgll_p;
  Idxs::HostMirror geo_c2n, cgll_c2n, cgll_io_c2n;
  mesh::make_cubedsphere(geo_p, geo_c2n, in.n);
  Int nerr = 0;
  mesh::make_cgll_from_geo(geo_p, geo_c2n, np, cgll_p, cgll_c2n);
  mesh::make_io_cgll_from_internal_cgll(cgll_p, cgll_c2n, cgll_io_c2n);
  { // Clip the mesh against itself and get the total area.
    const Real
      area = siqk::test::test_area_ot<geometry>(cgll_p, cgll_io_c2n,
                                                cgll_p, cgll_io_c2n),
      true_area = 4*M_PI,
      re = std::abs(area - true_area)/true_area;
    if (re >= 1e-10) {
      fprintf(stderr, "true area %1.4e mesh area %1.4e relerr %1.4e\n",
              true_area, area, re);
      ++nerr;
    }
  }
  {
    const Int ne = mesh::impl::check_elem_normal_against_sphere(
      cgll_p, cgll_io_c2n);
    if (ne) std::cerr << "FAIL: check_elem_normal_against_sphere\n";
    nerr += ne;
  }
  {
    IdxArray::HostMirror dglln2cglln;
    Idxs::HostMirror dgll_c2n;
    mesh::make_dgll_from_cgll(cgll_p, cgll_c2n, dglln2cglln, dgll_c2n);
    const Int np2 = szslice(cgll_c2n);
    Int ne = 0;
    for (Int ci = 0; ci < nslices(cgll_c2n); ++ci) {
      const auto cgll_cell = slice(cgll_c2n, ci);
      for (Int ni = 0; ni < siqk::square(np); ++ni)
        if (dglln2cglln[ci*np2 + ni] != cgll_cell[ni])
          ++ne;
    }
    if (ne) {
      nerr += ne;
      std::cerr << "FAIL: make_dgll_from_cgll\n";
    }
  }
  if ( ! in.fn_pre_out.empty()) {
    io::NetcdfWriter ncw(cgll_p, cgll_io_c2n, in.fn_pre_out + ".g");
    ncw.add_nodal_field("x");
    ncw.end_definition();
    const Int n = nslices(cgll_p);
    std::vector<Real> x(n), lat(n), lon(n);
    for (Int i = 0; i < n; ++i) {
      const auto p = slice(cgll_p, i);
      xyz2ll(p[0], p[1], p[2], lat[i], lon[i]);
    }
    ncw.advance_time_to(1);
    gallery::InitialCondition::init(
      gallery::InitialCondition::CosineBells,
      nslices(cgll_p), lat.data(), lon.data(), x.data());
    ncw.write_field("x", x.data());
    ncw.advance_time_to(1.5);
    gallery::InitialCondition::init(
      gallery::InitialCondition::SlottedCylinders,
      nslices(cgll_p), lat.data(), lon.data(), x.data());
    ncw.write_field("x", x.data());
    ncw.advance_time_to(2.5);
    gallery::InitialCondition::init(
      gallery::InitialCondition::CorrelatedCosineBells,
      nslices(cgll_p), lat.data(), lon.data(), x.data());
    ncw.write_field("x", x.data());
  }
  if (in.write_matlab) {
    write_matlab("cm", geo_p, geo_c2n);
    write_matlab("m", cgll_p, cgll_io_c2n);
    write_matlab("gll", cgll_p, cgll_c2n);
  }
  return nerr;
}

static Int test_time_int (const Input& in) {
  return timeint::test::test_ark( ! in.quiet);
}

Input::Input (Int argc, char** argv)
  : command(Command::test_make_cubedsphere), n(10), angle(M_PI*1e-1),
    write_matlab(false), quiet(false)
{
  for (Int i = 1; i < argc; ++i) {
    const std::string& token = argv[i];
    if (eq(token, "-c", "--command")) command = Command::from_string(argv[++i]);
    else if (eq(token, "-n")) n = atoi(argv[++i]);
    else if (eq(token, "-q", "--quiet")) quiet = true;
    else if (eq(token, "-m", "--write-matlab")) write_matlab = true;
    else if (eq(token, "-o", "--output-prefix")) fn_pre_out = argv[++i];
    else if (eq(token, "--angle")) angle = atof(argv[++i]);
  }

  if ( ! quiet) print(std::cout);
}

void Input::print (std::ostream& os) const {
  os << "command " << command << "\n"
     << "n (-n): " << n << "\n"
     << "write matlab (-m): " << write_matlab << "\n"
     << "angle (--angle): " << angle << "\n";
}

int main (int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  {
    Input in(argc, argv);
    Int nerr = 0;
    switch (in.command) {
    case Command::test_make_cubedsphere: nerr = test_make_cubedsphere(in); break;
    case Command::test_gll: nerr = test_gll(in); break;
    case Command::test_make_gll_mesh: nerr = test_make_gll_mesh(in); break;
    case Command::test_time_int: nerr = test_time_int(in); break;
    }
    std::cerr << (nerr ? "FAIL" : "PASS") << "ED\n";
  }
  Kokkos::finalize_all();
}
