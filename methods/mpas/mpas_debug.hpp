#ifndef INCLUDE_MPAS_DEBUG_HPP
#define INCLUDE_MPAS_DEBUG_HPP

#include <sys/time.h>

#include <cstdlib>
#include <iostream>
#include <sstream>

#include "mpas_types.hpp"

#ifdef MPAS_HAVE_MPI
# include <mpi.h>
#endif

namespace dbg {
static inline int getrank () {
  int pid = 0;
#ifdef MPAS_HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
#endif
  return pid;
}

static inline int getnrank () {
  int n = 0;
#ifdef MPAS_HAVE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &n);
#endif
  return n;
}

static void
write_array_impl (const Array1D<const int>* gids, const Array3D<const double>& d,
                  const int d0, const int d1, const int d2beg, const int d2end,
                  const std::string& filename) {
  const int rank = dbg::getrank(), nrank = dbg::getnrank();
  std::stringstream ss;
  ss << filename << "_" << nrank << "_" << rank << ".txt";
  FILE* fid = fopen(ss.str().c_str(), "wa");
  for (int d2 = d2beg; d2 < d2end; ++d2)
    fprintf(fid, "%d %23.16e\n", gids ? (*gids)(d2) : d2, d(d0,d1,d2));
  fclose(fid);
}

void write_array(const Array1D<const int>* gids, const Array3D<const double>& d,
                 const int d0, const int d1, const int d2beg, const int d2end,
                 const std::string& filename);

// Silence warning re: unused static functions.
inline void unused () {
  Array3D<const double> unused;
  write_array_impl(nullptr, unused, 0, 0, 0, 0, "");
}
} // namespace dbg

#ifndef pr
# ifdef MPAS_HAVE_MPI
#define pr(m) do {                                      \
    int _pid_ = dbg::getrank();                         \
    std::stringstream _ss_;                             \
    _ss_ << "pid " << _pid_ << " " << m << std::endl;   \
    std::cerr << _ss_.str();                            \
  } while (0)
#define pr0(m) do {                                     \
    int _pid_ = dbg::getrank();                         \
    MPI_Comm_rank(MPI_COMM_WORLD, &_pid_);              \
    if (_pid_ != 0) break;                              \
    std::stringstream _ss_;                             \
    _ss_ << "pid " << _pid_ << " " << m << std::endl;   \
    std::cerr << _ss_.str();                            \
  } while (0)
# else
#define pr(m) do {                                      \
    std::stringstream _ss_;                             \
    _ss_ << m << std::endl;                             \
    std::cerr << _ss_.str();                            \
  } while (0)
#define pr0(m) pr(m)
# endif
#define prc(m) pr(#m << " | " << (m))
#define pr0c(m) pr0(#m << " | " << (m))
#define puf(m) "(" << #m << " " << (m) << ")"
#define pu(m) << " " << puf(m)
#endif

template <typename T>
static void prarr (const std::string& name, const T* const v, const size_t n) {
  std::stringstream ss;
  ss << name << " = [";
  for (size_t i = 0; i < n; ++i) ss << " " << v[i];
  ss << "]";
  pr(ss.str());
}
template <typename Array>
static void prarr (const std::string& name, const Array& a) {
  prarr(name, a.data(), a.size());
}

#ifndef mprarr
# define mprarr(a) prarr(#a, a)
#endif

#endif
