#ifndef INCLUDE_MPAS_UTIL_HPP
#define INCLUDE_MPAS_UTIL_HPP

#include "mpas_types.hpp"
#include "mpas_siqp.hpp"

#ifdef _OPENMP
# include <omp.h>
#endif
#include <sys/time.h>

#define throw_if(condition, message) do {                               \
    if (condition) {                                                    \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": The condition:\n"       \
           << #condition "\nled to the exception\n" << message << "\n"; \
        throw std::logic_error(_ss_.str());                             \
    }                                                                   \
  } while (0)
#define require(condition) do {                                         \
    if ( ! (condition)) {                                               \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": FAIL:\n" << #condition  \
        << "\n";                                                        \
      throw std::logic_error(_ss_.str());                               \
    }                                                                   \
  } while (0)
#define require_msg(condition, message) do {                            \
    if ( ! (condition)) {                                               \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": FAIL:\n" << #condition  \
           << "\nmessage:\n" << message << "\n";                        \
      throw std::logic_error(_ss_.str());                               \
    }                                                                   \
  } while (0)

int mpas_omp_get_thread_num();
int mpas_omp_get_max_threads();

template <typename T> inline T square (const T& a) { return a*a; }
template <typename T> inline T cube (const T& a) { return a*a*a; }
template<typename T> inline T sign (const T& a) { return a >= 0 ? 1 : -1; }

struct FILECloser { void operator() (FILE* fh) { fclose(fh); } };
using FilePtr = std::unique_ptr<FILE, FILECloser>;

template <int nrow=3, typename T>
inline const T* slice (const T* v, const int idx) { return v + nrow*idx; }
template <int nrow=3, typename T>
inline T* slice (T* v, const int idx) { return v + nrow*idx; }

inline Real reldif (const Real& a, const Real& b) {
  return (a == 0 ?
          std::abs(b - a) :
          std::abs(b - a)/std::abs(a));
}

double urand();

bool eq(const std::string& a, const char* const b1, const char* const b2 = 0);

bool fileexists(const std::string& fname);

// C = alpha op(A) op(B) + beta C.
void dgemm(char transa, char transb, int m, int nrhs, int n, double alpha,
           const double* a, int lda, const double* b, int ldb, double beta,
           const double* c, int ldc);

// Cholesky factorization and solve.
int dpotrf(char uplo, int n, double* a, int lda);
int dpotrs(char uplo, int n, int nrhs, const double* a, int lda, double* bx,
           int ldb);
void dtrsm(char side, char uplo, char transa, char diag, int n, int nrhs,
           double alpha, const double* a, int lda, double* bx, int ldb);

// QR factorization functions.
int dtrtrs(char uplo, char trans, char diag, int n, int nrhs,
           double* a, int lda, double* b, int ldb);
// tau[min(m,n)], wrk[>= n]
int dgeqrf(int m, int n, double* a, int lda,
           double* tau, double* wrk, int iwrk);
int dorgqr(int m, int n, int k, double* a, int lda, double* tau,
           double* wrk, int iwrk);
// tau[min(m,n)], wrk[>= max(m,n)]
int dormqr(char side, char trans, int m, int n, int k, double* a, int lda,
           double* tau, double* c, int ldc, double* wrk, int iwrk);

bool lsinv(
  // Column-major matrix with n rows and ncol columns.
  Real* A, const int n, const int ncol,
  // Work arrays. tau, wrk are for BLAS calls; Q has size A.
  Real* tau, Real* wrk, Real* Q,
  // Inverse of operator, row-major.
  Real* B);

void lsinv_linear_shape_function(const int n, const Real* xs, const Real* ys, Real* B);

/* Solve min_x norm(A x - b, 2)^2
          st   c'x = d
   for mxn full-column-rank A, nxk x and b, nx1 c, 1xk d, m >= n. In the inputs
   and outputs, A, bx, dy are column major, and all are overwritten. bx contains
   b on input and x on output; dy, d and the Lagrange multipler.
     rwrk should have size at least
       m*n + m*k + std::min(m,n) + std::max(m,n).
*/
Int solve_ls_1eq(const int m, const int n, const int k, Real* A, Real* cs,
                 Real* bx, Real* dy, Real* rwrk);

template <int max_nhit_>
struct HitsAccumulator {
  static const int max_nhit = max_nhit_;
  int hits[max_nhit], nhit;

  HitsAccumulator () : nhit(0) {}

  bool operator() (const int i) {
    if (nhit == max_nhit) consolidate();
    assert(nhit <= max_nhit); // should not overflow buffer
    if (nhit == max_nhit) return false; // no room in buffer
    hits[nhit++] = i;
    return true;
  }

  void post (const Int* gids = nullptr) {
    if (gids)
      std::sort(hits, hits + nhit,
                [&] (const int i, const int j) { return gids[i] < gids[j]; });
    else
      std::sort(hits, hits + nhit);
    int prev = hits[0];
    for (int i = 1; i < nhit; ++i) {
      if (hits[i] == prev) hits[i] = -1;
      else prev = hits[i];
    }
  }

private:
  void consolidate () { /* No impl right now. */ }
};

struct Timer {
  enum Op { integrate, trajectory, reconstruct, search, clip, quadrature, scq,
            cedr, rest, total, ntimers };
  static const int max_nthreads = 64;

  static inline void init () {
#ifdef MPAS_TIME
    for (int i = 0; i < ntimers*max_nthreads; ++i) et[i] = 0;
#endif
  }

  static inline void start (const Op op) { start_impl(0, op); }
  static inline void stop (const Op op) { stop_impl(0, op); }

  static inline void start (const int tid, const Op op) {
#ifdef MPAS_TIME_INNER
    start_impl(tid, op);
#endif
  }
  static inline void stop (const int tid, const Op op) {
#ifdef MPAS_TIME_INNER
    stop_impl(tid, op);
#endif
  }

  static void print();

private:
#ifdef MPAS_TIME
  static timeval t_start[max_nthreads * ntimers];
  static double et[max_nthreads * ntimers];
#endif

  static inline void start_impl (const int tid, const Op op) {
#ifdef MPAS_TIME
    gettimeofday(&t_start[tid*ntimers + op], 0);
#endif
  }

  static inline void stop_impl (const int tid, const Op op) {
#ifdef MPAS_TIME
    timeval t2;
    gettimeofday(&t2, 0);
    const int idx = tid*ntimers + op;
    const timeval& t1 = t_start[idx];
    static const double us = 1.0e6;
    et[idx] += (t2.tv_sec*us + t2.tv_usec - t1.tv_sec*us - t1.tv_usec)/us;
#endif
  }
};

#endif
