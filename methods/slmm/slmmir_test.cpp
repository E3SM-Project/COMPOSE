#include <stdexcept>

#include "slmm_defs.hpp"
#include "slmm_mesh.hpp"
#include "slmm_debug.hpp"
using namespace slmm;

struct Command {
  enum Enum { test_make_cubedsphere };
  static Enum from_string (const std::string& s) {
    if (s == "test_make_cubedsphere") return test_make_cubedsphere;
    throw std::runtime_error(s + " is not a command.");
  }
};

struct Input {
  Command::Enum command;
  Int n;
  Real angle;
  bool write_matlab;

  Input(Int argc, char** argv);
  void print(std::ostream& os) const;
};

static Int test_make_cubedsphere (const Input& in) {
  AVec3s cp;
  AIdxs ce;
  mesh::make_cubedsphere(cp, ce, in.n);
  Int nerr = 0;
  {
    const Int ne = mesh::check_elem_normal_against_sphere(cp, ce);
    if (ne) std::cerr << "FAIL: check_elem_normal_against_sphere\n";
    nerr += ne;
  }
  if (in.write_matlab)
    write_matlab("cm", cp, ce);
  return nerr;
}

inline bool
eq (const std::string& a, const char* const b1, const char* const b2 = 0) {
  return (a == std::string(b1) || (b2 && a == std::string(b2)) ||
          a == std::string("-") + std::string(b1));
}

Input::Input (Int argc, char** argv)
  : command(Command::test_make_cubedsphere), n(10), angle(M_PI*1e-1),
    write_matlab(false)
{
  for (Int i = 1; i < argc; ++i) {
    const std::string& token = argv[i];
    if (eq(token, "-c", "--comand")) command = Command::from_string(argv[++i]);
    else if (eq(token, "-n")) n = atoi(argv[++i]);
    else if (eq(token, "-m", "--write-matlab")) write_matlab = true;
    else if (eq(token, "--angle")) angle = atof(argv[++i]);
  }

  print(std::cout);
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
    case Command::test_make_cubedsphere:
      nerr = test_make_cubedsphere(in);
      break;
    }
    std::cerr << (nerr ? "FAIL" : "PASS") << "ED\n";
  }
  Kokkos::finalize_all();
}
