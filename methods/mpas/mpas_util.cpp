#include "mpas_util.hpp"
#include "mpas_mesh.hpp"

double urand () { return rand() / ((double) RAND_MAX + 1.0); }

bool
eq (const std::string& a, const char* const b1, const char* const b2) {
  return (a == std::string(b1) || (b2 && a == std::string(b2)) ||
          a == std::string("-") + std::string(b1));
}

int mpas_omp_get_thread_num () {
#ifdef _OPENMP
  return omp_get_thread_num();
#endif
  return 0;
}

int mpas_omp_get_max_threads () {
#ifdef _OPENMP
  return omp_get_max_threads();
#endif
  return 1;  
}

bool fileexists (const std::string& fname) {
  auto* f = fopen(fname.c_str(), "r");
  if (f) {
    fclose(f);
    return true;
  }
  return false;
}

extern "C" {
  void dgemm_(const char* transa, const char* transb, const int* m,
              const int* n, const int* k, const double* alpha, const double* a,
              const int* lda, const double* b, const int* ldb,
              const double* beta, double* c, const int* ldc);
  void dpotrf_(const char* uplo, const int* n, double* a, const int* lda,
               int* info);
  void dpotrs_(const char* uplo, const int* n, const int* nrhs, const double* a,
               const int* lda, double* b, const int* ldb, int* info);
  void dtrsm_(const char* side, const char* uplo, const char* transa, const char* diag,
              const int* n, const int* nrhs, const double* alpha, const double* a,
              const int* lda, double* b, const int* ldb);
  void dtrtrs_(const char* uplo, const char* trans, const char* diag,
               const int* n, const int* nrhs, double* a, const int* lda,
               double* b, const int* ldb, int* info);
  void dgeqrf_(const int* m, const int* n, double* a, const int* lda,
               double* tau, double* wrk, int* iwrk, int* info);
  void dorgqr_(const int* m, const int* n, const int* k, double* a,
               const int* lda, double* tau, double* wrk, int* iwrk, int* info);
  void dormqr_(const char* side, const char* trans,
               const int* m, const int* n, const int* k,
               double* a, const int* lda,
               double* tau, double* c, const int* ldc,
               double* wrk, const int* iwrk, int* info);
  void dgetrf_(const int* m, const int* n, double* a, const int* lda, int* ipiv,
               int* info);
  void dgetrs_(const char* trans, const int* n, const int* nrhs, const double* a,
               const int* lda, const int* ipiv, double* b, const int* ldb, int* info);
}

// C = alpha op(A) op(B) + beta C
void dgemm (char transa, char transb, int m, int nrhs, int n, double alpha,
            const double* a, int lda, const double* b, int ldb, double beta,
            const double* c, int ldc) {
  dgemm_(&transa, &transb, &m, &nrhs, &n, &alpha, const_cast<double*>(a), &lda,
         const_cast<double*>(b), &ldb, &beta, const_cast<double*>(c), &ldc);
}

int dpotrf (char uplo, int n, double* a, int lda) {
  int info;
  dpotrf_(&uplo, &n, a, &lda, &info);
  return info;
}

int dpotrs (char uplo, int n, int nrhs, const double* a, int lda, double* bx,
            int ldb) {
  int info;
  dpotrs_(&uplo, &n, &nrhs, const_cast<double*>(a), &lda, bx, &ldb, &info);
  return info;
}

void dtrsm (char side, char uplo, char transa, char diag, int n, int nrhs,
            double alpha, const double* a, int lda, double* bx, int ldb) {
  dtrsm_(&side, &uplo, &transa, &diag, &n, &nrhs, &alpha,
         const_cast<double*>(a), &lda, bx, &ldb);
}

int dtrtrs (char uplo, char trans, char diag, int n, int nrhs,
            double* a, int lda, double* b, int ldb) {
  int info;
  dtrtrs_(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, &info);
  return info;
}

// tau[min(m,n)], wrk[>= n]
int dgeqrf (int m, int n, double* a, int lda,
            double* tau, double* wrk, int iwrk) {
  int info;
  dgeqrf_(&m, &n, a, &lda, tau, wrk, &iwrk, &info);
  return info;
}

int dorgqr (int m, int n, int k, double* a, int lda, double* tau,
            double* wrk, int iwrk) {
  int info;
  dorgqr_(&m, &n, &k, a, &lda, tau, wrk, &iwrk, &info);
  return info;
}

// tau[min(m,n)], wrk[>= max(m,n)]
int dormqr (char side, char trans, int m, int n, int k, double* a, int lda,
            double* tau, double* c, int ldc, double* wrk, int iwrk) {
  int info;
  dormqr_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, wrk, &iwrk, &info);
  return info;
}

int dgetrf (int m, int n, double* a, int lda, int* ipiv) {
  int info;
  dgetrf_(&m, &n, a, &lda, ipiv, &info);
  return info;
}

int dgetrs (char trans, int n, int nrhs, const double* a, int lda, const int* ipiv,
            double* bx, int ldb) {
  int info;
  dgetrs_(&trans, &n, &nrhs, a, &lda, ipiv, bx, &ldb, &info);
  return info;
}

bool lsinv (Real* A, const int n, const int ncol,
            Real* tau, Real* wrk, Real* Q, Real* B) {
  int info = dgeqrf(n, ncol, A, n, tau, wrk, n); assert(info == 0);
  std::copy(A, A + ncol*n, Q);
  info = dorgqr(n, ncol, ncol, Q, n, tau, wrk, n); assert(info == 0);
  Real* rhs = tau;
  bool ok = true;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < ncol; ++j) rhs[j] = Q[n*j+i];
    info = dtrtrs('U', 'N', 'N', ncol, 1, A, n, rhs, ncol);
    if (info) ok = false;
    for (int j = 0; j < ncol; ++j) B[n*j+i] = rhs[j];
  }
  return ok;
}

void lsinv_linear_shape_function (const int n, const Real* xs, const Real* ys, Real* B) {
  static const int ncol = 3;
  assert(n >= 3);
  Real A[3*max_1halo_slots];
  for (int i = 0; i < n; ++i) {
    A[n*0+i] = 1;
    A[n*1+i] = xs[i];
    A[n*2+i] = ys[i];
  }
  Real tau[ncol*max_1halo_slots], wrk[max_1halo_slots], Q[ncol*max_1halo_slots];
  lsinv(A, n, ncol, tau, wrk, Q, B);
}

Int solve_ls_1eq (const int m, const int n, const int k, Real* A, Real* cs,
                  Real* bx, Real* dy, Real* rwrk) {
  /* The Lagrangian is
         Lag = 1/2 x' A'A x - b' A x - mu (c' x - d),
     giving the system
         Lag_x  = A'A x - A' b - c mu = 0
         Lag_mu = c'x - d = 0,
     or
         [A'A -c] [x ] = [A'b]
         [c'    ] [mu]   [d  ].
     We have Q R = qr(A) and L = R'. Then
         [L          0] [L'  -inv(L) c     ] [x] = [L Q'b]
         [c'inv(L')  I] [0    c'inv(L L') c] [y]   [d    ].
     Hence
         1. [L          0] [a1] = [L Q'b]
            [c'inv(L')  I] [a2]   [d    ]
         => 1a. Solve L s = c
            1b. Set a1 = Q'b
            1c. Set a2 = d - s'a1.
         2. [L'  -inv(L) c     ] [x] = [a1]
            [0    c'inv(L L') c] [y]   [a2]
         => 2a. Set y = a2 / (s's)
            2b. Solve L' x = a1 + y s.
  */
  Real* const R = A; // alias
  Real* const Q = rwrk;
  Real* const tau = Q + m*n;
  Real* const a1 = tau + std::min(m,n);
  Real* const wrk = a1 + n*k;
  // Form Q R = qr(A).
  int info = dgeqrf(m, n, A, m, tau, wrk, std::max(m, n));
  if (info) return info;
  std::copy(A, A + m*n, Q);
  info = dorgqr(m, n, n, Q, m, tau, wrk, std::max(m, n));
  if (info) return info;
  // 1b
  dgemm('t', 'n', n, k, m, 1, Q, m, bx, m, 0, a1, n);
  // 1a
  info = dtrtrs('u', 't', 'n', n, 1, R, m, cs, n);
  if (info) return info;
  // 1c and 2a
  Real s2 = 0;
  for (int j = 0; j < n; ++j) s2 += square(cs[j]);
  for (int i = 0; i < k; ++i) {
    Real accum = 0;    
    for (int j = 0; j < n; ++j) accum += cs[j]*a1[n*i+j];
    dy[i] = (dy[i] - accum)/s2;
  }
  // 2b
  for (int i = 0; i < k; ++i)
    for (int j = 0; j < n; ++j)
      bx[m*i+j] = a1[n*i+j] + dy[i]*cs[j];
  info = dtrtrs('u', 'n', 'n', n, k, R, m, bx, m);
  return info;
}

#ifdef MPAS_TIME
timeval Timer::t_start[Timer::ntimers*Timer::max_nthreads];
double Timer::et[Timer::ntimers*Timer::max_nthreads];
#endif

# define tpr(op) do {                                             \
    printf("%-20s %10.3e %10.1f\n", #op, et[op], 100*et[op]/tot); \
  } while (0)
void Timer::print () {
#ifdef MPAS_TIME
    for (int j = 0; j < ntimers; ++j)
      for (int i = 1; i < max_nthreads; ++i)
        et[j] = std::max(et[j], et[i*ntimers + j]);
    const double tot = et[total];
    tpr(integrate);
      tpr(trajectory);
      tpr(reconstruct);
      tpr(scq);
#ifdef MPAS_TIME_INNER
        tpr(search); tpr(clip); tpr(quadrature);
#endif
      tpr(cedr);
    //tpr(rest);
    printf("%-20s %10.3e %10.1f\n", "total", et[total], 100.0);
#endif
  }
#undef tpr
