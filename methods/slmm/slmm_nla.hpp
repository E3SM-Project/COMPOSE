#ifndef INCLUDE_SLMM_NLA_HPP
#define INCLUDE_SLMM_NLA_HPP

#include <memory>

#include "slmm_array.hpp"
#include "slmm_util.hpp"

namespace slmm {

// C = alpha op(A) op(B) + beta C.
void dgemm(char transa, char transb, int m, int nrhs, int n, double alpha,
           const double* a, int lda, const double* b, int ldb, double beta,
           const double* c, int ldc);

// Cholesky factorization and solve.
int dpotrf(char uplo, int n, double* a, int lda);
int dpotrs(char uplo, int n, int nrhs, const double* a, int lda, double* bx,
           int ldb);
int dtrsm(char side, char uplo, char transa, char diag, int n, int nrhs,
          double alpha, const double* a, int lda, double* bx, int ldb);

// LU factorization and solve.
int dgetrf(int m, int n, double* a, int lda, int* ipiv);
int dgetrs(char trans, int n, int nrhs, const double* a, int lda,
           const int* ipiv, double* b, int ldb);

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

// SVD.
int dgesvd(char jobu, char jobvt, int m, int n, double* a, int lda,
           double* s, double* u, int ldu, double* vt, int ldvt,
           double* work, int lwork);

// y += A x for row-major, 4x4 A.
inline void matvec4(const double* __restrict__ a, const double* __restrict__  x,
                    double* __restrict__ y, const int nrhs) {
  assert(nrhs == 1);
  y[0] += a[ 0]*x[0] + a[ 1]*x[1] + a[ 2]*x[2] + a[ 3]*x[3];
  y[1] += a[ 4]*x[0] + a[ 5]*x[1] + a[ 6]*x[2] + a[ 7]*x[3];
  y[2] += a[ 8]*x[0] + a[ 9]*x[1] + a[10]*x[2] + a[11]*x[3];
  y[3] += a[12]*x[0] + a[13]*x[1] + a[14]*x[2] + a[15]*x[3];
}

// Row-major matrix.
void matvec(const Int m, const Int n, const Real* const a,
            const Real* const x, Real* const y);

// Row-major matrix.
void tmatvec(const Int m, const Int n, const Real* const a,
             const Real* const x, Real* const y);

// C = u C + v A B, for row-major A m by p and column-major B p by n, C m by n.
void matmult_rcc(const Int m, const Int n, const Int p,
                 const Real v, const Real* const A, const Real* const B,
                 const Real u, Real* const C);

// C =  u C + v A B, for column-major A, B, C.
void matmult_ccc(const Int m, const Int n, const Int p,
                 const Real v, const Real* const A, const Real* const B,
                 const Real u, Real* const C);

/* Solve
       [A  B] [X] = [C]
       [B' 0] [Y]   [D],
   where m >= n (= implies the constraints fully specify the solution) and
       A is m by m and s.p.d., B is m by n,
       C is m by k, D is n by k,
       X is m by k, Y is n by k.
   On input, A, B, CX, DY are column-major. On output, CX contains X, DY
   contains Y, and A is overwritten with work. rwrk has size m*n and iwrk has
   size n.
 */
Int solve_kkt(const Int m, const Int n, const Int k,
              Real* const A, const Real* const B,
              Real* const CX, Real* const DY,
              Real* const rwrk, Int* const iwrk);

Int test_solve_kkt();

// Form op = A \ B for spd A m by m, B m by n, both column-major. On output, A
// is overwritten with internal computations. B_op contains B on input and op on
// output.
Int form_ls_op(const Int m, const Int n, Real* const A, Real* const B_op);

Int test_form_ls_op();

// A is m by n and column-major.
void transpose(const Int m, const Int n, const Real* const A, Real* const At);

// All indices and sizes are relative to blocks except m() and n().
// Whether each block is row- or col-major is up to the caller.
// Each row's cols must be sorted.
template <typename ScalarT = Real, typename SizeT = int, typename IntT = int>
class BlockMatrix {
public:
  typedef ScalarT Scalar;
  typedef SizeT Size;
  typedef IntT Int;

  typedef BlockMatrix<Scalar, Size, Int> Me;

  // Don't need N, really, but it's handy for assertions/debugging.
  Int M_, N_, m_, n_;
  std::shared_ptr<Size> rowptr_p_;
  std::shared_ptr<Int> colidx_p_;
  std::shared_ptr<Scalar> d_p_;
  Size* rowptr_;
  Int* colidx_;
  Scalar* d_;

public:
  BlockMatrix ()
    : M_(0), m_(0), n_(0), rowptr_(nullptr), colidx_(nullptr), d_(nullptr)
  {}

  BlockMatrix (const Int M, const Int N, const Int m, const Int n,
               const Size* rowptr, const Int* colidx) {
    init(M, N, m, n, rowptr, colidx);
  }

  void init(const Int M, const Int N, const Int m, const Int n,
            const Size* rowptr, const Int* colidx);

  const Size* rowptr () const { return rowptr_; }
  const Int* colidx () const { return colidx_; }

  const Int M () const { return M_; }
  const Int N () const { return N_; }
  const Int m () const { return m_; }
  const Int n () const { return n_; }

  const Scalar* blockrow (const Int br) const {
    assert(br < M_);
    return d_ + rowptr_[br]*m_*n_;
  }
  Scalar* blockrow (const Int br) {
    return const_cast<Scalar*>(const_cast<const Me*>(this)->blockrow(br));
  }

  const Scalar* block (const Int br, const Int bc) const {
    assert(br < M_);
    assert(bc < N_);
    const Int* const beg = colidx_ + rowptr_[br];
    const Int* const end = colidx_ + rowptr_[br+1];
    const Int* const idx = std::lower_bound(beg, end, bc);
    if (idx == end) return nullptr;
    const Int i = static_cast<Int>(idx - colidx_);
    return d_ + i*m_*n_;
  }
  Scalar* block (const Int br, const Int bc) {
    return const_cast<Scalar*>(const_cast<const Me*>(this)->block(br, bc));
  }

  void zero () { for (Size i = 0; i < rowptr_[M_]*m_*n_; ++i) d_[i] = 0; }

  static void test();
};

struct QrFac {
  typedef std::shared_ptr<QrFac> Ptr;

  QrFac (const Int& m, const Int& n)
    : m_(m), n_(n), d_(2*m*n + std::min(m,n) + std::max(m,n)),
      A_(d_.data()), Q_(A_ + m*n), tau_(Q_ + m*n), wrk_(tau_ + std::min(m,n))
  {
    assert(m_ >= n_);
  }

  int m () const { return m_; }
  int n () const { return n_; }
  Real* A () { return A_; }

  int factor () {
    int info = dgeqrf(m_, n_, A_, m_, tau_, wrk_, std::max(m_, n_));
    if (info) return info;
    std::copy(A_, A_ + m_*n_, Q_);
    return dorgqr(m_, n_, n_, Q_, m_, tau_, wrk_, std::max(m_, n_));
  }

  int solve_R (Real* const c) const {
    return dtrtrs('u', 'n', 'n', n_, 1, A_, m_, c, m_);
  }

  int solve_Rt (Real* const c) const {
    return dtrtrs('u', 't', 'n', n_, 1, A_, m_, c, m_);
  }

  // c = Q'b
  int apply_Qt (const Real* const b, Real* const c, const Int ncol=1) const {
    dgemm('t', 'n', n_, ncol, m_, 1, Q_, m_, b, m_, 0, c, n_);
    return 0;
  }

  // c = Q b
  int apply_Q (const Real* const b, Real* const c, const Int ncol=1) const {
    dgemm('n', 'n', m_, ncol, n_, 1, Q_, m_, b, n_, 0, c, m_);
    return 0;
  }

  static int unittest(const int n);

private:
  Int m_, n_;
  Array<Real> d_;
  Real* const A_;
  Real* const Q_;
  Real* const tau_;
  Real* const wrk_;
};

} // namespace slmm

#endif
