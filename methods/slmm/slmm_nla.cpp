#include <limits>

#include "slmm_util.hpp"
#include "slmm_nla.hpp"

extern "C" {
  void dgemm_(const char* transa, const char* transb, const int* m,
              const int* n, const int* k, const double* alpha, const double* a,
              const int* lda, const double* b, const int* ldb,
              const double* beta, double* c, const int* ldc);
  void dpotrf_(const char* uplo, const int* n, double* a, const int* lda,
               int* info);
  void dpotrs_(const char* uplo, const int* n, const int* nrhs, const double* a,
               const int* lda, double* b, const int* ldb, int* info);
  void dgetrf_(const int* m, const int* n, double* a, const int* lda, int* ipiv,
               int* info);
  void dgetrs_(const char* trans, const int* n, const int* nrhs, const double* a,
               const int* lda, const int* ipiv, double* b, const int* ldb, int* info);
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
  void dgesvd_(char* jobu, char* jobvt, int* m, int* n, double* a, int* lda,
               double* s, double* u, int* ldu, double* vt, int* ldvt,
               double* work, int* lwork, int* info);
}

namespace slmm {

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
  dpotrs_(&uplo, &n, &nrhs, a, &lda, bx, &ldb, &info);
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

int dtrsm (char side, char uplo, char transa, char diag, int n, int nrhs,
            double alpha, const double* a, int lda, double* bx, int ldb) {
  dtrsm_(&side, &uplo, &transa, &diag, &n, &nrhs, &alpha,
         const_cast<double*>(a), &lda, bx, &ldb);
  return 0;
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

int dgesvd (char jobu, char jobvt, int m, int n, double* a, int lda,
            double* s, double* u, int ldu, double* vt, int ldvt,
            double* work, int lwork) {
  int info;
  dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
          work, &lwork, &info);
  return info;
}

// Row-major matrix.
void matvec (const Int m, const Int n, const Real* const a,
             const Real* const x, Real* const y) {
  for (Int i = 0; i < m; ++i) {
    Real accum = 0;
    const Real* const ai = a + i*n;
    for (Int j = 0; j < n; ++j) accum += ai[j]*x[j];
    y[i] = accum;
  }
}

// Row-major matrix.
void tmatvec (const Int m, const Int n, const Real* const a,
              const Real* const x, Real* const y) {
  for (Int i = 0; i < n; ++i) y[i] = 0;
  for (Int j = 0; j < m; ++j) {
    const Real* const aj = a + n*j;
    for (Int i = 0; i < n; ++i)
      y[i] += aj[i]*x[j];
  }
}

// C = u C + v A B, for row-major A m by p and column-major B p by n, C m by n.
void matmult_rcc (const Int m, const Int n, const Int p,
                  const Real v, const Real* const A, const Real* const B,
                  const Real u, Real* const C) {
  for (Int j = 0; j < n; ++j) {
    Real* const Ccol = C + m*j;
    const Real* const Bcol = B + p*j;
    if (u == 0)
      for (Int i = 0; i < m; ++i)
        Ccol[i] = 0;
    for (Int i = 0; i < m; ++i) {
      const Real* const Arow = A + p*i;
      Real accum = 0;
      for (Int k = 0; k < p; ++k)
        accum += Arow[k]*Bcol[k];
      Ccol[i] = u*Ccol[i] + v*accum;
    }
  }
}

// C =  u C + v A B, for column-major A, B, C.
void matmult_ccc (const Int m, const Int n, const Int p,
                  const Real v, const Real* const A, const Real* const B,
                  const Real u, Real* const C) {
  for (Int j = 0; j < n; ++j) {
    Real* const Ccol = C + m*j;
    if (u == 0)
      for (Int i = 0; i < m; ++i)
        Ccol[i] = 0;
    else
      for (Int i = 0; i < m; ++i)
        Ccol[i] *= u;
    const Real* const Bcol = B + p*j;
    for (Int k = 0; k < p; ++k) {
      const Real* const Acol = A + m*k;
      const Real bk = Bcol[k];
      for (Int i = 0; i < m; ++i)
        Ccol[i] += v*Acol[i]*bk;
    }
  }
}

/* Solve
       [A  B] [X] = [C]
       [B' 0] [Y]   [D],
   where m >= n (= implies the constraints fully specify the solution) and
       A is m by m, B is m by n,
       C is m by k, D is n by k,
       X is m by k, Y is n by k.
   On input, A, B, CX, DY are column-major. On output, CX contains X, DY
   contains Y, and A and B are overwritten with work. rwrk has size m*n and iwrk
   has size n.
 */
Int solve_kkt (const Int m, const Int n, const Int k,
               Real* const A, const Real* const B,
               Real* const CX, Real* const DY,
               Real* const rwrk, Int* const iwrk) {
  /* Algorithm. Rearranging,
        Y = (B'A"B)"(B'A"C - D)
          = (B'V)"(B'W - D)
        X = A"C - A"B (B'A"B)"(B'A"C - D)
          = W - V (B'V)"(B'W - D).
     Thus:
        factor L'L = A
        solve  L'L W = C          use CX in the code
        solve  L'L V = B          V is m by n; use V = rwrk
        factor P L U = B'V
        solve  P L U Y = D - B'W  use P = iwrk
        set    X = W + V Y.
  */
  assert(m >= n);
  if (dpotrf('L', m, A, m) != 0) return 1;
  if (dpotrs('L', m, k, A, m, CX, m) != 0) return 2;
  Real* const V = rwrk;
  std::copy(B, B + m*n, V);
  if (dpotrs('L', m, n, A, m, V, m) != 0) return 3;
  Real* const BV = A;
  matmult_rcc(n, n, m, 1, B, V, 0, BV);
  Int* const ipiv = iwrk;
  if (dgetrf(n, n, BV, n, ipiv) != 0) return 4;
  matmult_rcc(n, k, m, 1, B, CX, -1, DY);
  if (dgetrs('N', n, k, BV, n, ipiv, DY, n) != 0) return 5;
  matmult_ccc(m, k, n, -1, V, DY, 1, CX);
  return 0;
}

void transpose (const Int m, const Int n, const Real* const A, Real* const At) {
  for (Int i = 0; i < m; ++i)
    for (Int j = 0; j < n; ++j)
      At[n*i + j] = A[m*j + i];
}

// A = D'D
static void form_spd_A (const Int m, const Real* D, Real* const A) {
  matmult_rcc(m, m, m, 1, D, D, 0, A);
}

Int test_solve_kkt () {
  const Int sets[][3] = {{11,7,15}, {11,7,9}, {11,7,5}};
  Int nerr = 0;
  for (size_t iset = 0; iset < sizeof(sets)/sizeof(*sets); ++iset) {
    const Int m = sets[iset][0], n = sets[iset][1], k = sets[iset][2];
    Array<Real> Ac(m*m,0), Bc(m*n,0), X(m*k,0), Y(n*k,0), rwrk(m*m);
    Array<Int> iwrk(n);
    for (Int i = 0; i < m*m; ++i) rwrk[i] = urand() - 0.5;
    form_spd_A(m, rwrk.data(), Ac.data());
    for (Int i = 0; i < m; ++i)
      for (Int j = 0; j < m; ++j)
        if (Ac[i*m+j] != Ac[j*m+i]) ++nerr;  
    for (Int i = 0; i < m*n; ++i) Bc[i] = urand() - 0.5;
    for (Int i = 0; i < m*k; ++i) X[i] = urand() - 0.5;
    for (Int i = 0; i < n*k; ++i) Y[i] = urand() - 0.5;
    Array<Real> A(Ac), B(Bc), C(X), D(Y);
    nerr += solve_kkt(m, n, k, Ac.data(), Bc.data(), X.data(), Y.data(),
                      rwrk.data(), iwrk.data());
    Array<Real> Ct(m*k,0), Dt(n*k,0);
    matmult_ccc(m, k, m, 1, A.data(), X.data(), 0, Ct.data());
    matmult_ccc(m, k, n, 1, B.data(), Y.data(), 1, Ct.data());
    matmult_rcc(n, k, m, 1, B.data(), X.data(), 0, Dt.data());
    for (Int mat = 0; mat < 2; ++mat) {
      const Real* r, * t;
      Int sz;
      if (mat == 0) { r = C.data(); t = Ct.data(); sz = m*k; }
      else { r = D.data(); t = Dt.data(); sz = n*k; }
      Real num = 0, den = 0;
      for (Int i = 0; i < sz; ++i) num += square(t[i] - r[i]);
      for (Int i = 0; i < sz; ++i) den += square(r[i]);
      if (std::sqrt(num/den) > 1e4*std::numeric_limits<Real>::epsilon()) {
        pr(std::sqrt(num/den));
        ++nerr;
      }
    }
  }
  return nerr;
}

Int form_ls_op (const Int m, const Int n, Real* const A, Real* const B_op) {
  if (dpotrf('L', m, A, m) != 0) return 1;
  if (dpotrs('L', m, n, A, m, B_op, m) != 0) return 2;
  return 0;
}

Int test_form_ls_op () {
  const Int sets[][3] = {{11,7}, {11,20}};
  Int nerr = 0;
  for (size_t iset = 0; iset < sizeof(sets)/sizeof(*sets); ++iset) {
    const Int m = sets[iset][0], n = sets[iset][1];
    Array<Real> Ac(m*m,0), op(m*n,0), rwrk(m*m);
    for (Int i = 0; i < m*m; ++i) rwrk[i] = urand() - 0.5;
    form_spd_A(m, rwrk.data(), Ac.data());
    for (Int i = 0; i < m*n; ++i) op[i] = urand() - 0.5;
    Array<Real> A(Ac), B(op), rhs(m*n), tmp(m*n);
    nerr += form_ls_op(m, n, Ac.data(), op.data());
    transpose(m, n, op.data(), tmp.data()); // just to test this also
    transpose(n, m, tmp.data(), op.data());
    matmult_ccc(m, n, m, 1, A.data(), op.data(), 0, rhs.data());
    const Int sz = m*n;
    Real num = 0, den = 0;
    for (Int i = 0; i < sz; ++i) num += square(rhs[i] - B[i]);
    for (Int i = 0; i < sz; ++i) den += square(B[i]);
    if (std::sqrt(num/den) > 1e3*std::numeric_limits<Real>::epsilon()) {
      pr(std::sqrt(num/den));
      ++nerr;
    }
  }
  return nerr;
}

template <typename Scalar, typename Size, typename Int>
void BlockMatrix<Scalar, Size, Int>
::init (const Int M, const Int N, const Int m, const Int n,
        const Size* rowptr, const Int* colidx) {
  M_ = M; N_ = N; m_ = m; n_ = n;
  rowptr_p_ = std::shared_ptr<Size>(new Size[M_ + 1],
                                    std::default_delete<Size[]>());
  rowptr_ = rowptr_p_.get();
  memcpy(rowptr_, rowptr, (M_ + 1)*sizeof(Size));
  colidx_p_ = std::shared_ptr<Int>(new Int[rowptr[M_]],
                                   std::default_delete<Int[]>());
  colidx_ = colidx_p_.get();
  memcpy(colidx_, colidx, rowptr[M_]*sizeof(Int));
  d_p_ = std::shared_ptr<Scalar>(new Scalar[rowptr_[M_]*m_*n_],
                                 std::default_delete<Scalar[]>());
  d_ = d_p_.get();
}

template <typename Scalar, typename Size, typename Int>
void BlockMatrix<Scalar, Size, Int>::test () {
  static const Size rowptr[] = {0, 2, 3, 6, 8 };
  static const Size colidx[] = {0, 1, 1, 0, 2, 3, 1, 3};
  static const Int M = sizeof(rowptr)/sizeof(Size) - 1;
  static const int m = 3, n = 4;

  {
    BlockMatrix<Scalar, Size, Int> a(M, M, m, n, rowptr, colidx);

    assert(a.M() == M);
    assert(a.m() == m);
    assert(a.n() == n);

    const auto rowptr = a.rowptr();
    const auto colidx = a.colidx();
    for (Int r = 0, ctr = 1; r < a.M(); ++r) {
      Scalar* d = a.blockrow(r);
      for (Int j = 0; j < rowptr[r+1] - rowptr[r]; ++j, ++ctr) {
        for (Int i = 0; i < a.m()*a.n(); ++i)
          d[i] = ctr;
        d += a.m()*a.n();
      }
    }

    for (Int r = 0, ctr = 1; r < M; ++r)
      for (Int j = rowptr[r]; j < rowptr[r+1]; ++j, ++ctr) {
        Scalar const* const d = a.block(r, colidx[j]);
        assert(d);
        for (Int i = 0; i < m*n; ++i)
          assert(d[i] == ctr);
      }
  }
}

template class BlockMatrix<Real, Int, Int>;

int QrFac::unittest (const int n) {
  int nerr = 0;
  const auto m = n+1, mn = (n+1)*n;
  QrFac qr(m, n);
  auto A = qr.A();
  std::vector<Real> Acopy(mn), b(m), x(n), Ax_test(m);
  for (int i = 0; i < mn; ++i) Acopy[i] = A[i] = 2*urand() - 1;
  // Make the last row trivially dependent on the others to test m>1 but with
  // effectively an nxn problem.
  for (int j = 0; j < n; ++j) Acopy[m*j + n] = A[m*j + n] = A[m*j + n-1];
  qr.factor();
  for (int i = 0; i < n; ++i) b[i] = 2*urand() - 1;
  b[n] = b[n-1];
  qr.apply_Qt(b.data(), x.data());
  qr.apply_Q(x.data(), Ax_test.data());
  qr.solve_R(x.data());
  std::vector<Real> Ax(m);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < m; ++i)
      Ax[i] += Acopy[m*j + i]*x[j];
  Real num1 = 0, num2 = 0, den = 0;
  for (int i = 0; i < m; ++i) {
    auto d = Ax[i] - b[i];
    num1 += d*d;
    d = Ax_test[i] - b[i];
    num2 += d*d;
    den += b[i]*b[i];
  }
  auto re = std::sqrt(num1/den);
  if (re > 1e-13) ++nerr;
  re = std::sqrt(num2/den);
  if (re > 1e-13) ++nerr;
  return nerr;
}

} // namespace slmm
