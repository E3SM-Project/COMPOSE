#ifndef INCLUDE_CEDR_UTIL_HPP
#define INCLUDE_CEDR_UTIL_HPP

#include <sstream>

#include "cedr_kokkos.hpp"
#include "cedr_mpi.hpp"

namespace cedr {
namespace util {

template <typename T> KOKKOS_INLINE_FUNCTION constexpr
T square (const T& x) { return x*x; }

bool eq(const std::string& a, const char* const b1, const char* const b2 = 0);

// Uniform rand in [0, 1).
Real urand();

#define pr(m) do {                                      \
    int _pid_ = 0;                                      \
    MPI_Comm_rank(MPI_COMM_WORLD, &_pid_);              \
    std::stringstream _ss_;                             \
    _ss_.precision(15);                                 \
    _ss_ << "pid " << _pid_ << " " << m << std::endl;   \
    std::cerr << _ss_.str();                            \
  } while (0)
#define pr0(m) do {                                     \
    int _pid_; MPI_Comm_rank(MPI_COMM_WORLD, &_pid_);   \
    if (_pid_ != 0) break;                              \
    std::stringstream _ss_;                             \
    _ss_ << "pid " << _pid_ << " " << m << std::endl;   \
    std::cerr << _ss_.str();                            \
  } while (0)
#define prc(m) pr(#m << " | " << (m))
#define pr0c(m) pr0(#m << " | " << (m))
#define puf(m) "(" << #m << " " << (m) << ")"
#define pu(m) << " " << puf(m)
template <typename T>
void prarr (const std::string& name, const T* const v, const size_t n) {
  std::stringstream ss;
  ss.precision(15);
  ss << name << " = [";
  for (size_t i = 0; i < n; ++i) ss << " " << v[i];
  ss << "];";
  pr(ss.str());
}
#define mprarr(m) cedr::util::prarr(#m, m.data(), m.size())

#ifndef NDEBUG
# define cedr_assert(condition) do {                                    \
    if ( ! (condition)) {                                               \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": FAIL:\n" << #condition  \
        << "\n";                                                        \
        throw std::logic_error(_ss_.str());                             \
    }                                                                   \
  } while (0)
# define cedr_kernel_assert(condition) do {     \
    if ( ! (condition))                         \
      Kokkos::abort(#condition);                \
  } while (0)
#else
# define cedr_assert(condition)
# define cedr_kernel_assert(condition)
#endif
#define cedr_throw_if(condition, message) do {                          \
    if (condition) {                                                    \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": The condition:\n"       \
           << #condition "\nled to the exception\n" << message << "\n"; \
        throw std::logic_error(_ss_.str());                             \
    }                                                                   \
  } while (0)
#define cedr_kernel_throw_if(condition, message) do {                   \
    if (condition)                                                      \
      Kokkos::abort(#condition " led to the exception\n" message);      \
  } while (0)

inline Real reldif (const Real a, const Real b)
{ return std::abs(b - a)/std::max(std::abs(a), std::abs(b)); }

Real reldif(const Real* a, const Real* b, const Int n);

struct FILECloser { void operator() (FILE* fh) { fclose(fh); } };

}
}

#endif
