#include <stdexcept>

#include "slmmir.hpp"

struct Command {
  enum Enum { test_csl_density };
  static Enum from_string (const std::string& s) {
    if (s == "test_csl_density") return test_csl_density;
    throw std::runtime_error(s + " is not a command.");
  }
};

struct Input {
  Command::Enum command;

  Input(Int argc, char** argv);
  void print(std::ostream& os) const;
};

static Int test_csl_density (const Input& in) {
  Int nerr = 0;
  return nerr;
}

Input::Input (Int argc, char** argv)
  : command(Command::test_csl_density)
{
  for (Int i = 1; i < argc; ++i) {
    const std::string& token = argv[i];
    if (eq(token, "-c", "--comand")) command = Command::from_string(argv[++i]);
  }

  print(std::cout);
}

void Input::print (std::ostream& os) const {
  os << "command " << command << "\n"
     << "\n";
}

int main (int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  {
    Input in(argc, argv);
    Int nerr = 0;
    switch (in.command) {
    case Command::test_csl_density:
      nerr = test_csl_density(in);
      break;
    }
    std::cerr << (nerr ? "FAIL" : "PASS") << "ED\n";
  }
  Kokkos::finalize_all();
}
