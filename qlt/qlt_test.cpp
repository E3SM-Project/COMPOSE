#include "qlt.hpp"

#include <stdexcept>
#include <sstream>

#define throw_if(condition, message) do {                               \
    if (condition) {                                                    \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": The condition:\n"       \
           << #condition "\nled to the exception\n" << message << "\n"; \
        throw std::logic_error(_ss_.str());                             \
    }                                                                   \
  } while (0)

inline bool eq (const std::string& a, const char* const b1, const char* const b2 = 0) {
  return (a == std::string(b1) || (b2 && a == std::string(b2)) ||
          a == std::string("-") + std::string(b1));
}

struct InputParser {
  qlt::test::Input in;

  class ArgAdvancer {
    const int argc_;
    char const* const* argv_;
    int i_;
  public:
    ArgAdvancer (int argc, char** argv) : argc_(argc), argv_(argv), i_(1) {}
    const char* advance () {
      if (i_+1 >= argc_) throw_if(true, "Command line is missing an argument.");
      return argv_[++i_];
    }
    const char* token () const { return argv_[i_]; }
    void incr () { ++i_; }
    bool more () const { return i_ < argc_; }
  };
  
  InputParser (int argc, char** argv, const qlt::Parallel::Ptr& p) {
    in.unittest = false;
    in.write = false;
    in.ncells = 0;
    in.ntracers = 1;
    in.tracer_type = 0;
    in.nrepeat = 1;
    in.pseudorandom = false;
    in.verbose = false;
    for (ArgAdvancer aa(argc, argv); aa.more(); aa.incr()) {
      const char* token = aa.token();
      if (eq(token, "-t", "--unittest")) in.unittest = true;
      else if (eq(token, "-w", "--write")) in.write = true;
      else if (eq(token, "-nc", "--ncells")) in.ncells = std::atoi(aa.advance());
      else if (eq(token, "-nt", "--ntracers")) in.ntracers = std::atoi(aa.advance());
      else if (eq(token, "-tt", "--tracertype")) in.tracer_type = std::atoi(aa.advance());
      else if (eq(token, "-nr", "--nrepeat")) in.nrepeat = std::atoi(aa.advance());
      else if (eq(token, "--random")) in.pseudorandom = true;
      else if (eq(token, "-v", "--verbose")) in.verbose = true;
      else throw_if(true, "Invalid token " << token);
    }
    throw_if(in.tracer_type < 0 || in.tracer_type >= 4, "Tracer type is out of bounds [0, 3].");
    throw_if(in.ntracers < 1, "Number of tracers is < 1.");
  }

  void print (std::ostream& os) const {
    os << "ncells " << in.ncells
       << " nrepeat " << in.nrepeat;
    if (in.pseudorandom) os << " random";
    os << "\n";
  }
};

int main (int argc, char** argv) {
  int ret = 0;
  MPI_Init(&argc, &argv);
  auto p = qlt::make_parallel(MPI_COMM_WORLD);
  srand(p->rank());
  Kokkos::initialize(argc, argv);
  try {
    InputParser inp(argc, argv, p);
    if (p->amroot()) inp.print(std::cout);
    ret = qlt::test::run(p, inp.in);
    if (p->amroot()) std::cout << (ret != 0 ? "FAIL" : "PASS") << "\n";
  } catch (const std::exception& e) {
    if (p->amroot())
      std::cerr << e.what();
  }
  Kokkos::finalize_all();
  MPI_Finalize();
  return ret;
}
