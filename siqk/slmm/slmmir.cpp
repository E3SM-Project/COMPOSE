#include "slmm_defs.hpp"
#include "slmm_mesh.hpp"
#include "slmm_gll.hpp"
#include "slmm_io.hpp"
#include "slmm_time_int.hpp"
#include "slmm_gallery.hpp"
#include "slmm_debug.hpp"
using namespace slmm;

// -----------------------------------------------------------------------------
// NLA stuff taken from cflexp1 tr_gll. All of this needs to be rewritten for
// Kokkos. My plan is to get the program running end to end correctly, and then
// I'll go back and transition things to Kokkos and to running on the GPU.

template <typename T> class Array {
  T* p_;
  std::size_t n_, cap_;
public:
  Array () { init(); }
  Array(std::size_t n);
  Array(std::size_t n, const T& init);
  ~Array () { clear(); }
  // Initialize the object with the assumption that all variables are uninit'ed
  // prior to calling.
  void init();
  void clear();
  // optclear means optionally clear. The function has the semantics of
  // clearing, but it may not actually release the memory.
  void optclear_and_resize(std::size_t n);
  // _ft indicates first touch.
  void optclear_and_resize_ft(std::size_t n);
  void optclear_and_resize(std::size_t n, const T& i);
  void optclear_and_reserve(std::size_t n);
  void optclear_and_reserve_ft(std::size_t n);
  T& operator[] (std::size_t i) { return p_[i]; }
  const T& operator[] (std::size_t i) const { return p_[i]; }
  T& back () { return p_[n_-1]; }
  const T& back () const { return p_[n_-1]; }
  std::size_t size () const { return n_; }
  bool empty () const { return size() == 0; }
  T* data () const { return p_; }
  // This does not realloc; reserve must provide the necessary memory. It does
  // not throw, either. It asserts.
  void unsafe_push_back(const T& e);
  T* begin () { return p_; }
  T* end () { return p_ + n_; }
  void set (const T& v) { for (std::size_t i = 0; i < n_; ++i) p_[i] = v; }
};

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

  void init (const Int M, const Int N, const Int m, const Int n,
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

class FullMassMatrix {
  typedef BlockMatrix<Real, int, int> MT;

  int np_;
  MT m_;

public:
  typedef std::shared_ptr<FullMassMatrix> Ptr;
  
  FullMassMatrix () : np_(0) {}
  FullMassMatrix (const int nelem, const int np) { init(nelem, np); }

  void init(const int nelem, const int np);

  int np2 () const { return np_*np_; }
  int np4 () const { return np2()*np2(); }

  const Real* block(const int i) const;
  Real* block(const int i);

  const MT& get_M () const { return m_; }

  void factor();
  void solve(const int elem, Real* const bx, const int nrhs,
             const int ldbx) const;
};

class RemapData {
public:
  typedef std::shared_ptr<RemapData> Ptr;
  typedef BlockMatrix<Real, Int, Int> MT;
  typedef Array<Real> VT;
  typedef siqk::Octree<geometry, 10> Octree;

  // Full block-diag target-target mass matrix, factored.
  FullMassMatrix fmm_;
  // Search tree over Eulerian mesh.
  Octree ot_;
  // Target-source matrix.
  MT T_;
  // Jacobian(ref square -> sphere).
  RealArray::HostMirror Jt_;
  // Eulerian mesh basis function integrals.
  RealArray::HostMirror dgbfi_, cgbfi_;

public:
  // Set up.
  FullMassMatrix& fmm () { return fmm_; }
  Octree& octree () { return ot_; }
  MT& T () { return T_; }
  RealArray::HostMirror& Jt () { return Jt_; }
  RealArray::HostMirror& dgbfi () { return dgbfi_; }
  RealArray::HostMirror& cgbfi () { return cgbfi_; }

  // Apply.
  Int T_nrows () const { return T_.M()*T_.m(); }
  Int T_ncols () const { return T_.N()*T_.n(); }
  const Octree& octree () const { return ot_; }
  const ConstRealArray::HostMirror& Jt () const { return Jt_; }
  const ConstRealArray::HostMirror& dgbfi () const { return dgbfi_; }
  const ConstRealArray::HostMirror& cgbfi () const { return cgbfi_; }

  // y = T x.
  void apply_T(const Real* x, const int ldx, Real* y, const int ldy,
               const int nrhs) const;
  // y = T' x. Not needed in practice, but used in check().
  void apply_T_transp(const Real* x, const int ldx, Real* y, const int ldy,
                      const int nrhs) const;
  // x = M_full \ b.
  void solve_M_full(Real* bx, const int nrhs, const int ldxb) const;
  // y = R_full x
  void apply_R_full(const Real* x, const int ldx, Real* y, const int ldy,
                    const int nrhs) const;
  // y = R_lump x
  void apply_R_lump(const Real* x, const int ldx, Real* y, const int ldy,
                    const int nrhs) const;

  // Perform and print some checks. Each entry of these Jacobians is the
  // integral over the spherical quad of a basis function. So it's really more
  // than just a Jacobian.
  void check(const Real* Js, const Real* Jt) const;
  // If T is expected to be identical to M (analytically), check how close it
  // really is. Works only before 'factor' is called.
  void compare_MT() const;
};

template<typename T> inline void touch (T* const p, const size_t n,
                                        const T& init = T()) {
  // 1 KB should be a safe lower bound on page size. Touch enough to touch every
  // page; I don't think there's any need to touch more memory than that.
  for (size_t i = 0; i < n; i += 1024 / sizeof(T))
    p[i] = init;
  // Make sure the last part is touched.
  if (n) p[n-1] = init;
}
template<typename T> inline T*
allocn (const size_t n, const bool first_touch = false) {
  if ( ! n) return 0;
  T* p = new T[n];
  if (first_touch) touch(p, n);
  return p;
}
template<typename T> inline void deln (T*& p) {
  if (p) delete[] p;
  p = 0;
}
template<typename T> inline void deln_const (const T* p) {
  if (p) delete[] p;
}
template<typename T> inline void del (T*& p) {
  if (p) delete p;
  p = 0;
}

template<typename T>
inline void Array<T>::init () {
  n_ = cap_ = 0;
  p_ = 0;
}

template<typename T>
inline Array<T>::Array (std::size_t n)
  : p_(0), n_(0), cap_(0)
{ optclear_and_resize(n); }

template<typename T>
inline Array<T>::Array (std::size_t n, const T& init)
  : p_(0), n_(0), cap_(0)
{ optclear_and_resize(n, init); }

template<typename T>
inline void Array<T>::clear () {
  n_ = cap_ = 0;
  deln(p_);
}

template<typename T>
inline void Array<T>::optclear_and_reserve (std::size_t n) {
  n_ = 0;
  if (n <= cap_) return;
  clear();
  p_ = allocn<T>(n);
  cap_ = n;
}

template<typename T>
inline void Array<T>::optclear_and_reserve_ft (std::size_t n) {
  n_ = 0;
  if (n <= cap_) return;
  clear();
  p_ = allocn<T>(n, true);
  cap_ = n;
}

template<typename T>
inline void Array<T>::optclear_and_resize (std::size_t n) {
  if (n <= cap_) {
    n_ = n;
    return;
  }
  optclear_and_reserve(n);
  n_ = n;
}

template<typename T>
inline void Array<T>::optclear_and_resize_ft (std::size_t n) {
  if (n <= cap_) {
    n_ = n;
    return;
  }
  optclear_and_reserve_ft(n);
  n_ = n;
}

template<typename T>
inline void Array<T>::optclear_and_resize (std::size_t n, const T& init) {
  optclear_and_resize(n);
  for (std::size_t i = 0; i < n_; ++i)
    memcpy(p_ + i, &init, sizeof(init));
}

template<typename T>
inline void Array<T>::unsafe_push_back (const T& e) {
  assert(n_ < cap_);
  p_[n_++] = e;
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

extern "C" {
  void dgemm_(const char* transa, const char* transb, const int* m,
              const int* n, const int* k, const double* alpha, const double* a,
              const int* lda, const double* b, const int* ldb,
              const double* beta, double* c, const int* ldc);
  void dpotrf_(const char* uplo, const int* n, double* a, const int* lda,
               int* info);
  void dpotrs_(const char* uplo, const int* n, const int* nrhs, const double* a,
               const int* lda, double* b, const int* ldb, int* info);
}

inline void dgemm (
  char transa, char transb, int m, int nrhs, int n, double alpha,
  const double* a, int lda, const double* b, int ldb, double beta,
  const double* c, int ldc)
{
  dgemm_(&transa, &transb, &m, &nrhs, &n, &alpha, const_cast<double*>(a), &lda,
         const_cast<double*>(b), &ldb, &beta, const_cast<double*>(c), &ldc);
}

void FullMassMatrix::init (const int nelem, const int np) {
  np_ = np;
  Array<int> rowptr(nelem + 1), colidx(nelem);
  for (int i = 0; i < nelem; ++i) {
    rowptr[i] = i;
    colidx[i] = i;
  }
  rowptr[nelem] = nelem;
  m_.init(nelem, nelem, np2(), np2(), rowptr.data(), colidx.data());
  m_.zero();
  assert(m_.m() == np2() && m_.n() == np2());
  assert(m_.M() == m_.N() && m_.M() == nelem);
  assert(m_.blockrow(0) + np4() == m_.blockrow(1));
}

const double* FullMassMatrix::block (const int i) const {
  assert(m_.blockrow(i) - m_.blockrow(0) == i*np4());
  return m_.blockrow(i);
}
double* FullMassMatrix::block (const int i) {
  return const_cast<double*>(const_cast<const MT&>(m_).blockrow(i));
}

void FullMassMatrix::factor () {
  const int n = np2();
# pragma omp parallel for
  for (int i = 0; i < m_.M(); ++i) {
    double* const d = block(i);
    const char uplo = 'L';
    int info;
    dpotrf_(&uplo, &n, d, &n, &info);
    if (info != 0) {
      fprintf(stderr, "M() %d i %d info %d\n", m_.M(), i, info);
      fprintf(stderr, "a = [");
      for (int c = 0; c < n; ++c) {
        for (int r = 0; r < n; ++r)
          fprintf(stderr, " %1.15e", d[n*c + r]);
        fprintf(stderr, ";");
      }
      fprintf(stderr, "];\n");
    }
    assert(info == 0);
  }
}

void FullMassMatrix::
solve (const int elem, double* const bx, const int nrhs, const int ldbx) const {
  const int n = np2();
  const double* const d = block(elem);
  const char uplo = 'L';
  int info;
  dpotrs_(&uplo, &n, &nrhs, const_cast<double*>(d), &n, bx, &ldbx, &info);
  assert(info == 0);
}

void RemapData::apply_T (const double* x, const int ldx, double* y,
                         const int ldy, const int nrhs) const {
  const MT::Scalar* const d = T_.blockrow(0);
  const MT::Size* const rowptr = T_.rowptr();
  const MT::Int* const colidx = T_.colidx();
# pragma omp parallel
  {
    const MT::Int n = T_.N()*T_.n();
#   pragma omp for
    for (MT::Int i = 0; i < n; ++i)
      y[i] = 0;
#   pragma omp for
    for (MT::Size br = 0; br < T_.M(); ++br)
      for (MT::Int j = rowptr[br]; j < rowptr[br+1]; ++j) {
        const MT::Int bc = colidx[j];
        const MT::Scalar* const b = d + j*T_.m()*T_.n();
        dgemm('t', 'n', T_.m(), nrhs, T_.n(), 1, b, T_.m(), x + bc*T_.n(), ldx,
              1, y + br*T_.m(), ldy);
      }
  }
}

void RemapData::apply_T_transp (const double* x, const int ldx, double* y,
                                const int ldy, const int nrhs) const {
  const MT::Scalar* const d = T_.blockrow(0);
  const MT::Size* const rowptr = T_.rowptr();
  const MT::Int* const colidx = T_.colidx();
  for (MT::Int i = 0, n = T_.M()*T_.m(); i < n; ++i)
    y[i] = 0;
  for (MT::Size br = 0; br < T_.M(); ++br)
    for (MT::Int j = rowptr[br]; j < rowptr[br+1]; ++j) {
      const MT::Int bc = colidx[j];
      const MT::Scalar* const b = d + j*T_.m()*T_.n();
      dgemm('n', 'n', T_.m(), nrhs, T_.n(), 1, b, T_.m(), x + br*T_.m(), ldx, 1,
            y + bc*T_.n(), ldy);
    }
}

void RemapData::solve_M_full (double* bx, const int nrhs,
                              const int ldxb) const {
# pragma omp parallel for
  for (MT::Int br = 0; br < T_.M(); ++br)
    fmm_.solve(br, bx + br*fmm_.np2(), nrhs, ldxb);
}

void RemapData::apply_R_full (const double* x, const int ldx, double* y,
                              const int ldy, const int nrhs) const {
  const MT::Int n = T_nrows();
  apply_T(x, n, y, n, 1);
  solve_M_full(y, 1, n);
}

static void report (const std::string label, const Real* const x_t,
                    const Real* const x, const Int n) {
  Real me = 0, den = 0;
  for (Int i = 0; i < n; ++i) {
    me = std::max(me, std::abs(x[i] - x_t[i]));
    den = std::max(den, std::abs(x_t[i]));
  }
  printf("> RemapData %21s: %1.3e\n", label.c_str(), me/den);
}

void RemapData::check (const Real* Js, const Real* Jt) const {
  const int n = T_nrows();
  // This routine assumes T is nxn.
  Array<double> e(n), x(n), y(n);
  e.set(1);

  memcpy(x.data(), Jt, n*sizeof(Real));
  solve_M_full(x.data(), 1, n);
  report("M_full \\ Jt = e", e.data(), x.data(), n);

  apply_T_transp(e.data(), n, x.data(), n, 1);
  report("e' T = Js'", Js, x.data(), n);

  apply_T(e.data(), n, x.data(), n, 1);
  report("T e = Jt", Jt, x.data(), n);

  apply_R_full(e.data(), n, x.data(), n, 1);
  report("[ct]     R_full e = e", e.data(), x.data(), n);

  memcpy(x.data(), Jt, n*sizeof(Real));
  solve_M_full(x.data(), 1, n);
  apply_T_transp(x.data(), n, y.data(), n, 1);
  report("[cv] Jt' R_full = Js'", Js, y.data(), n);
}

void RemapData::compare_MT () const {
  Real diag_num = 0, diag_den = 0;
  const auto& M = fmm_.get_M();
  const auto& T = T_;
  assert(M.M() == T.M());
  assert(M.m() == T.m());
  for (Int br = 0; br < T.M(); ++br) {
    const auto Mb = M.block(br, br);
    const auto Tb = T.block(br, br);
    for (Int k = 0; k < square(M.m()); ++k) {
      diag_num += square(Tb[k] - Mb[k]);
      diag_den += square(Mb[k]);
    }
  }
  printf("> rd(M,T) %1.3e\n", std::sqrt(diag_num/diag_den));
}

// -----------------------------------------------------------------------------
//   fwd = forward: The mesh at t_{n-1} is the departure mesh and is integrated
// forward in time. It is the source mesh.
//   bwd = backward: The mesh at t_n is the departure mesh and is integrated
// backward in time. It is the target mesh.
//   R = M \ T. M is the mass matrix. T is the mixed mass matrix mapping source
// to target.

// Some debug and code stuff.
namespace {
class Debug {
  int index_;
  std::string filename_;
  bool on_;

public:
  Debug ()
    : index_(1), filename_("dbgout.m"), on_(true)
  {
#ifdef SLMM_DEBUG
    FILE* fid = fopen(filename_.c_str(), "w");
    fclose(fid);
#endif
  }

  void advance () { ++index_; }

  void set_on (const bool set) { on_ = set; }

  template <typename CV3s>
  void write_p (const std::string& name, const CV3s& p) {
#ifdef SLMM_DEBUG
    if ( ! on_) return;
    FILE* fid = fopen(filename_.c_str(), "a");
    fprintf(fid, "%s{%d} = [", name.c_str(), index_);
    for (Int ip = 0; ip < nslices(p); ++ip)
      fprintf(fid, " %1.15e %1.15e %1.15e;", p(ip,0), p(ip,1), p(ip,2));
    fprintf(fid, "].';\n");
    fclose(fid);
#endif
  }

  template <typename CIs>
  void write_c2n (const std::string& name, const CIs& e) {
#ifdef SLMM_DEBUG
    if ( ! on_) return;
    FILE* fid = fopen(filename_.c_str(), "a");
    fprintf(fid, "%s{%d} = [", name.c_str(), index_);
    for (Int ie = 0; ie < nslices(e); ++ie) {
      for (Int k = 0; k < szslice(e); ++k)
        fprintf(fid, " %d", e(ie,k)+1);
      fprintf(fid, ";");
    }
    fprintf(fid, "].';\n");
    fclose(fid);
#endif
  }

  void write (const std::string& name, const BlockMatrix<Real, Size, Int>& m) {
#ifdef SLMM_DEBUG
    if ( ! on_) return;
    FILE* fid = fopen(filename_.c_str(), "a");
    fprintf(fid, "tmp = [");
    const Size* rowptr = m.rowptr();
    const Int* colidx = m.colidx();
    for (Int R = 0; R < m.M(); ++R)
      for (Int J = rowptr[R]; J < rowptr[R+1]; ++J) {
        const Int C = colidx[J];
        const Real* const block = m.block(R, C);
        for (Int r = 0, k = 0; r < m.m(); ++r)
          for (Int c = 0; c < m.n(); ++c, ++k)
            fprintf(fid, "%d %d %1.15e\n", m.m()*R + r + 1,
                    m.n()*C + c + 1, block[k]);
      }
    fprintf(fid, "];\n");
    fprintf(fid, "%s{%d} = sparse(tmp(:,1),tmp(:,2),tmp(:,3),%d,%d);\n",
            name.c_str(), index_, m.M()*m.m(), m.N()*m.n());
    fclose(fid);
#endif
  }

  void write (const std::string& name, const Real* const a, const Int n) {
#ifdef SLMM_DEBUG
    if ( ! on_) return;
    FILE* fid = fopen(filename_.c_str(), "a");
    fprintf(fid, "%s{%d} = [", name.c_str(), index_);
    for (Int i = 0; i < n; ++i)
      fprintf(fid, " %1.15e", a[i]);
    fprintf(fid, "].';\n");
    fclose(fid);
#endif
  }
};
static Debug gdbg;

class Timer {
public:
  enum Op { ts_setup, ts, ts_integrate, ts_remap, ts_rest, ts_error,
            ts_remap_T, ts_remap_node_jac,
            ts_remap_T_geometry, ts_remap_T_crs, ts_remap_T_fill,
            total, NTIMERS };
  static inline void init () {
#ifdef SLMM_TIME
    for (int i = 0; i < NTIMERS; ++i) et_[i] = 0;
#endif
  }
  static inline void start (const Op op) {
#ifdef SLMM_TIME
    gettimeofday(&t_start_[op], 0);
#endif
  }
  static inline void stop (const Op op) {
#ifdef SLMM_TIME
    timeval t2;
    gettimeofday(&t2, 0);
    const timeval& t1 = t_start_[op];
    static const double us = 1.0e6;
    et_[op] += (t2.tv_sec*us + t2.tv_usec - t1.tv_sec*us - t1.tv_usec)/us;
#endif
  }
# define tpr(op) do {                                                   \
    printf("%-20s %10.3e %10.1f\n", #op, et_[op], 100*et_[op]/tot);      \
  } while (0)
  static void print () {
#ifdef SLMM_TIME
    const double tot = et_[total];
    tpr(ts_setup); tpr(ts); tpr(ts_integrate); tpr(ts_remap);
    tpr(ts_remap_T); tpr(ts_remap_T_geometry); tpr(ts_remap_T_crs);
    tpr(ts_remap_T_fill); tpr(ts_remap_node_jac); tpr(ts_rest);
    tpr(ts_error);
    printf("%-20s %10.3e %10.1f\n", "total", et_[total], 100.0);
#endif
  }
#undef tpr
private:
#ifdef SLMM_TIME
  static timeval t_start_[NTIMERS];
  static double et_[NTIMERS];
#endif
};
#ifdef SLMM_TIME
timeval Timer::t_start_[Timer::NTIMERS];
double Timer::et_[Timer::NTIMERS];
#endif
} // anon namespace

static constexpr Int max_nvert = 8;
static constexpr Int max_hits = 25; // Covers at least a 2-halo.

class MeshIntegrator {
protected:
  std::vector<Real> ll_;
public:
  MeshIntegrator (const Int nnodes)
    : ll_(2*nnodes)
  {}
  virtual ~MeshIntegrator () {}
  std::vector<Real>& get_ll () { return ll_; }
  // Must be called from inside ||{}.
  virtual void integrate(const Real ts, const Real tf, Vec3s::HostMirror& p) =0;
};

template<typename OdeFn>
class MeshIntegratorWithOdeFn : public MeshIntegrator {
  std::vector<timeint::Workspace> ws_;
  std::vector<Real> initial_step_;
  bool use_xyz_form_;

public:
  MeshIntegratorWithOdeFn (const Int nnodes, const bool use_xyz_form = false)
    : MeshIntegrator(nnodes), initial_step_(nnodes, 1e-3),
      use_xyz_form_(use_xyz_form)
  {}

  virtual void integrate (const Real ts, const Real tf, Vec3s::HostMirror& p) {
    const Int nn = nslices(p);
    assert(2*nn == static_cast<Int>(ll_.size()));
    ws_.resize(omp_get_max_threads());
#   pragma omp parallel for schedule(static, 4)
    for (Int i = 0; i < nn; ++i) {
      const int tid = omp_get_thread_num();

      // Our primary interest in these numerical experiments is order of
      // accuracy when the flow field is exact. Hence here we use extremely
      // tight error tolerances.
      timeint::Options opts;
      opts.set_abs_tol(std::numeric_limits<Real>::epsilon());
      opts.set_rel_tol(1e2*std::numeric_limits<Real>::epsilon());
      opts.set_initial_step(initial_step_[i]);

      timeint::Info info;
      OdeFn fun;
      fun.set_xyz_form(use_xyz_form_);
      if ( ! use_xyz_form_) {
        Real lli[] = {ll_[2*i], ll_[2*i+1]};
        timeint::ark45(opts, fun, lli, 2, ts, tf, ws_[tid], &info);
        auto n = slice(p, i);
        ll2xyz(lli[0], lli[1], n[0], n[1], n[2]);
      } else {
        Real u[3];
        ll2xyz(ll_[2*i], ll_[2*i+1], u[0], u[1], u[2]);
        timeint::ark45(opts, fun, u, 3, ts, tf, ws_[tid], &info);
        geometry::normalize(u);
        auto n = slice(p, i);
        for (Int j = 0; j < 3; ++j) n[j] = u[j];
      }
      initial_step_[i] = info.good_initial_step;
    }
  }
};

class MeshRotator : public MeshIntegrator {
  Vec3s::HostMirror p_;
  Real axis_[3];

public:
  MeshRotator (const ConstVec3s::HostMirror& p)
    : MeshIntegrator(nslices(p))
  {
    axis_[0] = 0.2; axis_[1] = 0.7; axis_[2] = 1;
    geometry::normalize(axis_);
    ko::resize(p_, nslices(p), szslice(p));
    ko::deep_copy(p_, p);
  }

  virtual void integrate (const Real ts, const Real tf, Vec3s::HostMirror& p) {
    const Int nn = nslices(p);
    assert(2*nn == static_cast<Int>(ll_.size()));
    const Real
      period = day2sec(12),
      a = 2*M_PI*(tf - ts)/period;
    Real r[9];
    form_rotation(axis_, a, r);
#   pragma omp parallel for
    for (Int i = 0; i < nn; ++i) {
      auto n = slice(p_, i);
      const Real x = n[0], y = n[1], z = n[2];
      n = slice(p, i);
      n[0] = r[0]*x + r[1]*y + r[2]*z;
      n[1] = r[3]*x + r[4]*y + r[5]*z;
      n[2] = r[6]*x + r[7]*y + r[8]*z;
    }
  }
};

struct MeshIntegratorFactory : public gallery::WindFieldType {
  static std::shared_ptr<MeshIntegrator>
  create (const std::string& ode, const bool use_xyz_form,
          const ConstVec3s::HostMirror& p)
  { return create(from_string(ode), use_xyz_form, p); }

  static std::shared_ptr<MeshIntegrator>
  create (const Enum& ode, const bool use_xyz_form,
          const ConstVec3s::HostMirror& p) {
    const Int nnodes = nslices(p);
    switch (ode) {
    case Dcmip1d3ll:
      return std::make_shared<MeshIntegratorWithOdeFn<
        gallery::Dcmip1d3llOdeFn> >(nnodes, use_xyz_form);
    case NonDivergentWindField:
      return std::make_shared<MeshIntegratorWithOdeFn<
        gallery::NonDivergentWindField> >(nnodes, use_xyz_form);
    case DivergentWindField:
      return std::make_shared<MeshIntegratorWithOdeFn<
        gallery::DivergentWindField> >(nnodes, use_xyz_form);
    case NonDivergentWindFieldHack:
      return std::make_shared<MeshIntegratorWithOdeFn<
        gallery::NonDivergentWindFieldHack> >(nnodes, use_xyz_form);
    case Rotate:
      return std::make_shared<MeshRotator>(p);
    default:
      assert(0);
    }
  }
};

struct IntegrateOptions {
  enum Enum { fwd, bwd, test_looa };
  Enum stepping;
  bool d2c; // Each step, and in error, convert dgll <-> cgll.
};

struct Input {
  std::string output_fn, ode, initial_condition, program_name;
  Real T;
  Int ne, nsteps, write_every, monotone_type, np, tq_order;
  bool debug, write_matlab;
  bool xyz_form; // Integrate in (x,y,z) space instead of (lat,lon).
  IntegrateOptions integrate_options;

  Input(Int argc, char** argv);
  void print(std::ostream& os) const;
};

// _s is start and _e is end.
struct Output {
  Real
    l2_err, max_err, mass_s, mass_e, min_s, max_s, min_e, max_e,
    et_timestep,
    mass_gll_s, mass_gll_e;
};

struct RemapOptions {
  Int np, monotone_type;

  RemapOptions ()
    : np(4), monotone_type(0)
  {}
};

struct Mesh {
  Int np, tq_order;
  Vec3s::HostMirror geo_p, geo_nml, cgll_p;
  Idxs::HostMirror geo_c2n, geo_c2nml, cgll_c2n, dgll_c2n, cgll_io_c2n;
  IdxArray::HostMirror dglln2cglln;
};

static void copy_vertices (
  const ConstVec3s::HostMirror& p, const ConstIdxs::HostMirror& c2n,
  const Int ci, Real* ps)
{
  const auto cell = slice(c2n, ci);
  for (Int i = 0; i < szslice(c2n); ++i) {
    const auto n = slice(p, cell[i]);
    for (Int k = 0; k < 3; ++k) ps[k] = n[k];
    ps += 3;
  }
}

static void calc_node_jacobians (
  const Mesh& m, const ConstVec3s::HostMirror& p, RealArray::HostMirror& J_dg)
{
  const Int np2 = square(m.np);
  ko::resize(J_dg, nslices(m.geo_c2n)*np2);
  GLL gll;
  const Real* gll_x, * gll_wt;
  gll.get_coef(m.np, gll_x, gll_wt);
# pragma omp parallel for
  for (Int ci = 0; ci < nslices(m.geo_c2n); ++ci) {
    const auto cell = slice(m.geo_c2n, ci);
    for (Int j = 0, basis_idx = 0; j < m.np; ++j) {
      const Real b = 0.5*(gll_x[j] + 1);
      for (Int i = 0; i < m.np; ++i, ++basis_idx) {
        const Real a = 0.5*(gll_x[i] + 1);
        Real J[9];
        siqk::sqr::impl::calc_Jacobian(p, cell, a, b, J);
        geometry::cross(J, J+3, J+6);
        const Real jac = std::sqrt(geometry::norm2(J+6));
        J_dg(ci*np2 + basis_idx) = jac;
      }
    }
  }
}

static void calc_basis_function_integrals (
  const Int np, const Int tq_order, const ConstVec3s::HostMirror& p,
  const ConstIdxs::HostMirror& c2n, RealArray::HostMirror& dgbfi)
{
  const Int np2 = square(np);
  ko::resize(dgbfi, nslices(c2n)*np2);
  ko::deep_copy(dgbfi, 0);
  siqk::TriangleQuadrature tq;
  siqk::RawConstVec3s tq_bary;
  siqk::RawConstArray tq_w;
  tq.get_coef(tq_order, tq_bary, tq_w);
  const Int nq = len(tq_w);
  GLL gll;
# pragma omp parallel for
  for (Int ci = 0; ci < nslices(c2n); ++ci) { // cell
    Real ps[12];
    copy_vertices(p, c2n, ci, ps);
    const auto cell = slice(c2n, ci);
    for (Int k = 1; k <= 2; ++k) // 2 triangles per quad cell
      for (Int q = 0; q < nq; ++q) { // quad point
        Real sphere_coord[3];
        const Real jac = geometry::calc_tri_jacobian(
          ps, ps+3*k, ps+3*(k+1), slice(tq_bary, q), sphere_coord);
        Real gj[GLL::max_np], gi[GLL::max_np]; {
          Real a, b;
          siqk::sqr::calc_sphere_to_ref(p, cell, sphere_coord, a, b);
          gll.eval(np, b, gj);
          gll.eval(np, a, gi);
        }
        const Real d0 = 0.5 * tq_w[q] * jac;
        for (Int j = 0, basis_idx = 0; j < np; ++j) { // along ref y dir
          const Real d1 = d0 * gj[j];
          for (Int i = 0; i < np; ++i, ++basis_idx) // along ref x dir
            dgbfi(ci*np2 + basis_idx) += d1 * gi[i];
        }
      }
  }
}

static void calc_basis_function_integrals (
  const Mesh& m, const ConstVec3s::HostMirror& p, RealArray::HostMirror& dgbfi,
  RealArray::HostMirror& cgbfi)
{
  calc_basis_function_integrals(m.np, m.tq_order, p, m.geo_c2n, dgbfi);
  ko::resize(cgbfi, nslices(m.cgll_p));
  ko::deep_copy(cgbfi, 0);
  for (Int i = 0; i < len(m.dglln2cglln); ++i)
    cgbfi(m.dglln2cglln(i)) += dgbfi(i);
}

static void calc_gll_basis_function_integrals (
  const Mesh& m, const ConstVec3s::HostMirror& p, RealArray::HostMirror& J_dg)
{
  const Int np2 = square(m.np);
  ko::resize(J_dg, nslices(m.geo_c2n)*np2);
  GLL gll;
  const Real* gll_x, * gll_wt;
  gll.get_coef(m.np, gll_x, gll_wt);
# pragma omp parallel for
  for (Int ci = 0; ci < nslices(m.geo_c2n); ++ci) {
    const auto cell = slice(m.geo_c2n, ci);
    for (Int j = 0, basis_idx = 0; j < m.np; ++j) {
      const Real b = 0.5*(gll_x[j] + 1);
      for (Int i = 0; i < m.np; ++i, ++basis_idx) {
        const Real a = 0.5*(gll_x[i] + 1);
        Real J[9];
        siqk::sqr::impl::calc_Jacobian(p, cell, a, b, J);
        geometry::cross(J, J+3, J+6);
        const Real jac = std::sqrt(geometry::norm2(J+6));
        // Product of weights is the integral of the 2D basis function on the
        // ref square. Multiply by Jacobian of the map bilinear quad ->
        // sphere. Since this is GLL quadrature, there's exactly one quadrature
        // point.
        J_dg(ci*np2 + basis_idx) = 0.25 * jac * gll_wt[i] * gll_wt[j];
      }
    }
  }
}

static void map_cgll2dgll (
  const IdxArray::HostMirror& dglln2cglln, const Real* const cg_data,
  Real* const dg_data)
{
# pragma omp parallel for
  for (Int i = 0; i < len(dglln2cglln); ++i)
    dg_data[i] = cg_data[dglln2cglln[i]];
}

static void map_dgll2cgll (
  const IdxArray::HostMirror& dglln2cglln,
  const ConstRealArray::HostMirror& dgbfi,
  const ConstRealArray::HostMirror& cgbfi,
  const Real* const dg_data, Real* const cg_data, const Int cnn)
{
  for (Int i = 0; i < cnn; ++i) cg_data[i] = 0;
  for (Int i = 0; i < len(dglln2cglln); ++i) {
    const Int i_cgll = dglln2cglln(i);
    cg_data[i_cgll] += (dgbfi(i) / cgbfi(i_cgll)) * dg_data[i];
  }
}

static void calc_M_fwd (const Mesh& m, RemapData& rd) {
  const auto& p = m.geo_p;
  const auto& c2n = m.geo_c2n;
  auto& fmm = rd.fmm();
  fmm.init(nslices(c2n), m.np);
  const Int np = m.np, np2 = square(np);
  siqk::TriangleQuadrature tq;
  siqk::RawConstVec3s tq_bary;
  siqk::RawConstArray tq_w;
  tq.get_coef(m.tq_order, tq_bary, tq_w);
  const Int nq = len(tq_w);
  GLL gll;
# pragma omp parallel for
  for (Int ci = 0; ci < nslices(c2n); ++ci) {
    Real ps[12];
    copy_vertices(p, c2n, ci, ps);
    const auto cell = slice(c2n, ci);
    Real* block = fmm.block(ci);
    for (Int k = 1; k <= 2; ++k)
      for (Int q = 0; q < nq; ++q) {
        Real sphere_coord[3];
        const Real jac = geometry::calc_tri_jacobian(
          ps, ps+3*k, ps+3*(k+1), slice(tq_bary, q), sphere_coord);
        Real gj[GLL::max_np], gi[GLL::max_np]; {
          Real a, b;
          siqk::sqr::calc_sphere_to_ref(p, cell, sphere_coord, a, b);
          gll.eval(np, b, gj);
          gll.eval(np, a, gi);
        }
        const Real d0 = 0.5 * tq_w[q] * jac;
        for (Int aj = 0, a_basis_idx = 0; aj < np; ++aj) {
          const Real d1 = d0 * gj[aj];
          for (Int ai = 0; ai < np; ++ai, ++a_basis_idx) {
            const Real d2 = d1 * gi[ai];
            for (Int bj = 0, b_basis_idx = 0; bj < np; ++bj) {
              const Real d3 = d2 * gj[bj];
              for (Int bi = 0; bi < np; ++bi, ++b_basis_idx) {
                if (b_basis_idx < a_basis_idx) continue;
                const Real d = d3 * gi[bi];
                block[np2*a_basis_idx + b_basis_idx] += d;
                if (a_basis_idx != b_basis_idx)
                  block[np2*b_basis_idx + a_basis_idx] += d;
              }
            }
          }
        }
      }
  }
  gdbg.write("M", rd.fmm().get_M());
  //fmm.factor();
}

class CountIntersectionsFunctor {
protected:
  const siqk::sh::Mesh<ko::HostSpace>& cm_;
  const ConstVec3s::HostMirror p_;
  const ConstIdxs::HostMirror e_;
  Int hits_[max_hits];
  Int k_, nh_;

public:
  CountIntersectionsFunctor (
    const siqk::sh::Mesh<ko::HostSpace>& cm, const ConstVec3s::HostMirror& p,
    const ConstIdxs::HostMirror& c2n)
    : cm_(cm), p_(p), e_(c2n), nh_(0)
  {}

  void reset (const Int clipped_ci) {
    k_ = clipped_ci;
    nh_ = 0;
  }

  void operator() (const Int clip_ci) {
    // Check whether we've clipped against this polygon before and there was a
    // non-0 intersection.
    for (Int i = 0; i < nh_; ++i)
      if (hits_[i] == clip_ci)
        return;
    // We have not, so do the intersection.
    Int no = 0;
    {
      // Area of all overlapping regions.
      // In and out vertex lists.
      Real buf[9*max_nvert];
      siqk::RawVec3s
        vi(buf, max_nvert, 3),
        vo(buf + 3*max_nvert, max_nvert, 3),
        wrk(buf + 6*max_nvert, max_nvert, 3);
      Int ni;
      ni = 0;
      for (Int i = 0; i < szslice(e_); ++i) {
        if (e_(k_,i) == -1) break;
        geometry::copy(slice(vi, i), slice(p_, e_(k_,i)));
        ++ni;
      }
      siqk::sh::clip_against_poly<geometry>(cm_, clip_ci, vi, ni, vo, no, wrk);
    }
    if (no) {
      // Non-0 intersection, so record.
      if (nh_ == max_hits) Kokkos::abort("max_hits is too small.");
      hits_[nh_++] = clip_ci;
    }
  }

  Int get_nhits () const { return nh_; }
  const Int* get_hits () const { return hits_; }
};

static void calc_T_pattern_fwd (
  const Mesh& m, const ConstVec3s::HostMirror& depart_p,
  const RemapData::Octree& ot, std::vector<Int>& rowptr,
  std::vector<Int>& colidx)
{
  Timer::start(Timer::ts_remap_T_geometry);
  const Int ncell = nslices(m.geo_c2n);
  Idxs::HostMirror hits("hits", ncell, max_hits);
  {
    siqk::sh::Mesh<ko::HostSpace> cm;
    cm.p = m.geo_p; cm.e = m.geo_c2n; cm.nml = m.geo_nml; cm.en = m.geo_c2nml;
#   pragma omp parallel for schedule(static, 20)
    for (Int ci = 0; ci < ncell; ++ci) {
      Real bb[6];
      RemapData::Octree::calc_bb(depart_p, slice(m.geo_c2n, ci),
                                 szslice(m.geo_c2n), bb);
      CountIntersectionsFunctor cif(cm, depart_p, m.geo_c2n);
      cif.reset(ci);
      ot.apply(bb, cif);
      const Int* ci_hits = cif.get_hits();
      const Int hin = cif.get_nhits();
      for (Int hi = 0; hi < hin; ++hi)
        hits(ci, hi) = ci_hits[hi];
      if (hin < max_hits)
        hits(ci, hin) = -1;
    }
  }
  Timer::stop(Timer::ts_remap_T_geometry); Timer::start(Timer::ts_remap_T_crs);
  // Need to form transpose of the matrix that is most naturally created by the
  // above if using CRS format.
  rowptr.resize(ncell + 1, 0);
  for (Int ci = 0; ci < ncell; ++ci)
    for (Int hi = 0; hi < max_hits; ++hi) {
      if (hits(ci, hi) == -1) break;
      ++rowptr[hits(ci, hi) + 1];
    }
  // Cumsum.
  for (Int ci = 1; ci <= ncell; ++ci)
    rowptr[ci] += rowptr[ci-1];
  colidx.resize(rowptr[ncell]);
  // Shift up 1.
  for (Int ci = ncell; ci > 0; --ci)
    rowptr[ci] = rowptr[ci-1];
  for (Int ci = 0; ci < ncell; ++ci)
    for (Int hi = 0; hi < max_hits; ++hi) {
      const Int row = hits(ci, hi);
      if (row == -1) break;
      colidx[rowptr[row+1]] = ci;
      ++rowptr[row+1];
    }
# pragma omp parallel for
  for (Int ci = 0; ci < ncell; ++ci)
    std::sort(colidx.data() + rowptr[ci], colidx.data() + rowptr[ci+1]);
  Timer::stop(Timer::ts_remap_T_crs);
}

static void fill_T_fwd (const Mesh& m, const ConstVec3s::HostMirror& depart_p,
                        RemapData::MT& T) {
  const Int ncell = nslices(m.geo_c2n);
  const Int np = m.np, np2 = square(np), np4 = square(np2);
  siqk::TriangleQuadrature tq;
  siqk::RawConstVec3s tq_bary;
  siqk::RawConstArray tq_w;
  tq.get_coef(m.tq_order, tq_bary, tq_w);
  const Int nq = len(tq_w);
  GLL gll;
  const Size* rowptr = T.rowptr();
  const Int* colidx = T.colidx();
  siqk::sh::Mesh<ko::HostSpace> cm;
  cm.p = m.geo_p; cm.e = m.geo_c2n; cm.nml = m.geo_nml; cm.en = m.geo_c2nml;
# pragma omp parallel for schedule(static, 1)
  for (Int tci = 0; tci < ncell; ++tci) {
    Real* block = T.blockrow(tci);
    const auto tcell = slice(m.geo_c2n, tci);
    for (Int cj = rowptr[tci]; cj < rowptr[tci+1]; ++cj) {
      const Int sci = colidx[cj];
      const auto scell = slice(m.geo_c2n, sci);
      Real buf[9*max_nvert];
      siqk::RawVec3s
        vi(buf, max_nvert, 3),
        vo(buf + 3*max_nvert, max_nvert, 3),
        wrk(buf + 6*max_nvert, max_nvert, 3);
      Int ni = 0, no;
      for (Int i = 0; i < szslice(m.geo_c2n); ++i) {
        if (scell[i] == -1) break;
        geometry::copy(slice(vi, i), slice(depart_p, scell[i]));
        ++ni;
      }
      siqk::sh::clip_against_poly<geometry>(cm, tci, vi, ni, vo, no, wrk);
      assert(no);
      {
        for (Int i = 0; i < np4; ++i) block[i] = 0;
        for (Int ktri = 1; ktri < no-1; ++ktri) // triangles in vo
          for (Int q = 0; q < nq; ++q) { // quad point
            Real sphere_coord[3];
            const Real jac = geometry::calc_tri_jacobian(
              slice(vo,0), slice(vo,ktri), slice(vo,ktri+1), slice(tq_bary, q),
              sphere_coord);
            Real tgj[GLL::max_np], tgi[GLL::max_np],
              sgj[GLL::max_np], sgi[GLL::max_np];
            {
              Real ta, tb, sa, sb;
              siqk::sqr::calc_sphere_to_ref(m.geo_p, tcell, sphere_coord,
                                            ta, tb);
              siqk::sqr::calc_sphere_to_ref(depart_p, scell, sphere_coord,
                                            sa, sb);
              gll.eval(np, tb, tgj);
              gll.eval(np, ta, tgi);
              gll.eval(np, sb, sgj);
              gll.eval(np, sa, sgi);
            }
            const Real d0 = 0.5 * tq_w[q] * jac;
            for (Int tj = 0, t_basis_idx = 0; tj < np; ++tj) {
              const Real d1 = d0 * tgj[tj];
              for (Int ti = 0; ti < np; ++ti, ++t_basis_idx) {
                const Real d2 = d1 * tgi[ti];
                for (Int sj = 0, s_basis_idx = 0; sj < np; ++sj) {
                  const Real d3 = d2 * sgj[sj];
                  for (Int si = 0; si < np; ++si, ++s_basis_idx) {
                    const Real d = d3 * sgi[si];
                    block[np2*t_basis_idx + s_basis_idx] += d;
                  }
                }
              }
            }
          }
      }
      block += np4;
    }
  }
}

static void calc_T_fwd (const Mesh& m, const Vec3s::HostMirror& depart_p,
                        RemapData& rd)
{
  { // Build T's sparse matrix nonzero pattern.
    std::vector<Int> rowptr, colidx;
    calc_T_pattern_fwd(m, depart_p, rd.octree(), rowptr, colidx);
    const Int N = len(rowptr)-1, n = square(m.np);
    rd.T().init(N, N, n, n, rowptr.data(), colidx.data());
  }
  Timer::start(Timer::ts_remap_T_fill);
  fill_T_fwd(m, depart_p, rd.T());
  Timer::stop(Timer::ts_remap_T_fill);
}

// On input, src_tracer is rho*tracer. On output, it is just the updated
// tracer. Density is removed for output and error checking.
static void remap (
  RemapData& rd, const Mesh& m, const Vec3s::HostMirror& depart_p,
  Real* const src_tracer, Real* const tgt_tracer, const Int ntracers,
  Real* const src_density, Real* const tgt_density,
  // If in_dgll, we're working in DGLL space the whole time; otherwise, we're
  // doing CGLL -> DGLL -> remap -> CGLL. If in_dgll, wrk can be null;
  // otherwise, it must have length >= 2 dnn.
  const bool in_dgll, Real* const wrk)
{
  // For debugging and analysis, factor here.
  static bool first = true; if (first) {
    //rd.compare_MT();
    rd.fmm().factor();
    first = false;
  }

  Timer::start(Timer::ts_remap_T);
  const Int dnn = len(m.dglln2cglln), cnn = nslices(m.cgll_p),
    len = in_dgll ? dnn : cnn;
  calc_T_fwd(m, depart_p, rd);
  Timer::stop(Timer::ts_remap_T); Timer::start(Timer::ts_remap_node_jac);
  RealArray::HostMirror Js;
  calc_node_jacobians(m, depart_p, Js);
  Timer::stop(Timer::ts_remap_node_jac);

  for (Int ti = 0; ti < ntracers; ++ti) {
    Real* src, * tgt;
    if (in_dgll) {
      src = src_tracer + ti*len;
      tgt = tgt_tracer + ti*len;
    } else {
      src = wrk;
      tgt = wrk + dnn;
      map_cgll2dgll(m.dglln2cglln, src_tracer + ti*len, src);
    }
    // Adjust density according to the flow. At this point, the tracer field has
    // density in it.
#   pragma omp parallel for
    for (Int i = 0; i < dnn; ++i) {
      const Real q = rd.Jt()[i]/Js[i];
      src[i] *= q;
    }
    // L2 project.
    rd.apply_R_full(src, dnn, tgt, dnn, 1);
    if ( ! in_dgll)
      map_dgll2cgll(m.dglln2cglln, rd.dgbfi(), rd.cgbfi(), tgt,
                    tgt_tracer + ti*len, cnn);
  }

  {
    Real* src, * tgt;
    if (in_dgll) {
      src = src_density;
      tgt = tgt_density;
    } else {
      src = wrk;
      tgt = wrk + dnn;
      map_cgll2dgll(m.dglln2cglln, src_density, src);
    }
#   pragma omp parallel for
    for (Int i = 0; i < dnn; ++i) {
      const Real q = rd.Jt()[i]/Js[i];
      src[i] *= q;
    }
    rd.apply_R_full(src, dnn, tgt, dnn, 1);
    if ( ! in_dgll)
      map_dgll2cgll(m.dglln2cglln, rd.dgbfi(), rd.cgbfi(), tgt,
                    tgt_density, cnn);
  }

  // For output, remove density from tracer field.
  for (Int ti = 0; ti < ntracers; ++ti) {
#   pragma omp parallel for
    for (Int i = 0; i < len; ++i)
      tgt_tracer[ti*len + i] /= tgt_density[i];
  }
}

static void print_error (
  const Mesh& m, const ConstRealArray::HostMirror& J_gll, const bool in_dgll,
  const Real* const fs, const Real* const ds,
  const Real* const fe, const Real* const de, Output& out)
{
  Real l2_num = 0, l2_den = 0, max_num = 0, max_den = 0;
  out.max_s = -1e300; out.min_s = 1e300;
  out.max_e = -1e300; out.min_e = 1e300;
  out.mass_s = 0; out.mass_e = 0;
  out.mass_gll_s = 0; out.mass_gll_e = 0;
  siqk::TriangleQuadrature tq;
  siqk::RawConstVec3s tq_bary;
  siqk::RawConstArray tq_w;
  tq.get_coef(m.tq_order, tq_bary, tq_w);
  const Int nq = len(tq_w);
  GLL gll;
  const auto& c2n = m.geo_c2n;
  const auto& p = m.geo_p;
  const Int np = m.np, np2 = square(np);
  // GLL mass conservation.
  for (Int ci = 0; ci < nslices(c2n); ++ci)
    for (Int j = 0, basis_idx = 0; j < m.np; ++j)
      for (Int i = 0; i < m.np; ++i, ++basis_idx) {
        const Int k = ci*np2 + basis_idx;
        const Real w = J_gll[k];
        const Int idx = in_dgll ? k : m.cgll_c2n(ci, basis_idx);
        out.mass_gll_s += w * ds[idx];
        out.mass_gll_e += w * de[idx];
      }
  // Mass conservation wrt quadrature approximation of exact integrals.
  for (Int ci = 0; ci < nslices(c2n); ++ci) {
    const auto cell = slice(c2n, ci);
    Real ps[12];
    copy_vertices(p, c2n, ci, ps);
    for (Int k = 1; k <= 2; ++k)
      for (Int q = 0; q < nq; ++q) {
        Real sphere_coord[3];
        const Real jac = geometry::calc_tri_jacobian(
          ps, ps+3*k, ps+3*(k+1), slice(tq_bary, q), sphere_coord);
        Real gj[GLL::max_np], gi[GLL::max_np]; {
          Real a, b;
          siqk::sqr::calc_sphere_to_ref(p, cell, sphere_coord, a, b);
          gll.eval(np, b, gj);
          gll.eval(np, a, gi);
        }
        const Real d0 = 0.5 * tq_w[q] * jac;
        for (Int j = 0, basis_idx = 0; j < np; ++j) {
          const Real d1 = d0 * gj[j];
          for (Int i = 0; i < np; ++i, ++basis_idx) {
            const Int k = ci*np2 + basis_idx;
            const Int idx = in_dgll ? k : m.cgll_c2n(ci, basis_idx);
            const Real w = d1 * gi[i];
            const Real e = fe[idx] - fs[idx];
            out.mass_s += w * ds[idx];
            out.mass_e += w * de[idx];
            l2_num += w * square(e);
            l2_den += w * square(fs[idx]);
            max_num = std::max(max_num, std::abs(e));
            max_den = std::max(max_den, std::abs(fs[idx]));
            out.min_s = std::min(out.min_s, fs[idx]);
            out.max_s = std::max(out.max_s, fs[idx]);
            out.min_e = std::min(out.min_e, fe[idx]);
            out.max_e = std::max(out.max_e, fe[idx]);
          }
        }
      }
  }
  out.l2_err = std::sqrt(l2_num/l2_den);
  out.max_err = max_num/max_den;
  printf("> re l2 %9.3e max %9.3e\n", out.l2_err, out.max_err);
  printf("> [cv] re %10.3e\n", reldif(out.mass_s, out.mass_e));
  printf("> [cv gll] re %10.3e\n", reldif(out.mass_gll_s, out.mass_gll_e));
  printf("> [mo] min %10.3e %10.3e [%10.3e] max %10.3e %10.3e [%10.3e]\n",
         out.min_s, out.min_e, out.min_e - out.min_s,
         out.max_s, out.max_e, out.max_e - out.max_s);
}

static void print_one_liner (const Input& in, const Output& out) {
  std::cout << "<OL> method " << in.integrate_options.stepping
            << " ode " << in.ode << " ic " << in.initial_condition
            << " T " << in.T << " np " << in.np << " ne " << in.ne
            << " tq " << in.tq_order << " nsteps " << in.nsteps
            << " mono " << in.monotone_type;
  printf(" re l2 %9.3e max %9.3e", out.l2_err, out.max_err);
  printf(" cv re %9.3e", reldif(out.mass_s, out.mass_e));
  printf(" cvgll re %9.3e", reldif(out.mass_gll_s, out.mass_gll_e));
  printf(" mo min %9.3e %9.3e %9.3e max %9.3e %9.3e %9.3e",
         out.min_s, out.min_e, out.min_e - out.min_s,
         out.max_s, out.max_e, out.max_e - out.max_s);
  printf(" et ts %9.3e nthr %d", out.et_timestep, omp_get_max_threads());
  std::cout << " prog " << in.program_name;
  std::cout << " xyz " << in.xyz_form;
  std::cout << " d2c " << in.integrate_options.d2c;
  std::cout << "\n";
}

static void init_mesh (const Int np, const Int tq_order, const Int ne,
                       Mesh& m) {
  m.np = np;
  m.tq_order = tq_order;
  mesh::make_cubedsphere(m.geo_p, m.geo_c2n, ne);
  mesh::make_cgll_from_geo(m.geo_p, m.geo_c2n, np, m.cgll_p, m.cgll_c2n);
  mesh::make_dgll_from_cgll(m.cgll_p, m.cgll_c2n, m.dglln2cglln, m.dgll_c2n);
  mesh::make_io_cgll_from_internal_cgll(m.cgll_p, m.cgll_c2n, m.cgll_io_c2n);
  {
    siqk::sh::Mesh<ko::HostSpace> sm; sm.p = m.geo_p; sm.e = m.geo_c2n;
    siqk::test::fill_normals<geometry>(sm);
    m.geo_nml = sm.nml; m.geo_c2nml = sm.en;
  }
}

//   A bit of complication in this routine is opts.d2c. The natural thing to do
// is to work in DGLL space the whole time, except possibly when writing to the
// netcdf file. However, we need to mimic the intended application: at each
// step, we get CGLL fields, convert to DGLL, L2 project, then convert back. If
// opts.d2c, mimic this behavior; if ! opts.d2c, stay in DGLL space the whole
// time except when writing to the file. We support both behaviors so we can
// analyze the impact of going back and forth on accuracy.
static void integrate (
  const Mesh& m, const std::shared_ptr<MeshIntegrator>& mi,
  const RemapOptions& ro, const Real T, const Int nsteps,
  gallery::InitialCondition::Shape ic, const std::string& out_fn,
  const Int write_every, const IntegrateOptions opts, Output& out)
{
  Timer::start(Timer::ts_setup);
  const Int dnn = len(m.dglln2cglln), cnn = nslices(m.cgll_p),
    len = opts.d2c ? cnn : dnn;

  // Initialize I/O.
  std::shared_ptr<io::NetcdfWriter> ncw;
  if (write_every > 0) {
    ncw = std::make_shared<io::NetcdfWriter>(
      m.cgll_p, m.cgll_io_c2n, out_fn + ".g", ro.np, ro.monotone_type);
    ncw->add_nodal_field("tracer");
    ncw->add_nodal_field("density");
    ncw->end_definition();
  }

  // Eulerian mesh remap data.
  RealArray::HostMirror Jt_gll;
  calc_gll_basis_function_integrals(m, m.geo_p, Jt_gll);
  RemapData rd;
  calc_M_fwd(m, rd);
  rd.octree().init(m.geo_p, m.geo_c2n);
  calc_node_jacobians(m, m.cgll_p, rd.Jt());
  calc_basis_function_integrals(m, m.geo_p, rd.dgbfi(), rd.cgbfi());

  // Initialize data and workspace.
  std::vector<Real> tracer[2], density[2];
  std::vector<Real>* tracer_p[2], * density_p[2];
  for (Int i = 0; i < 2; ++i) {
    tracer[i].resize(len);
    tracer_p[i] = &tracer[i];
    density[i].resize(len);
    density_p[i] = &density[i];
  }
  for (Int k = 0; k < 2; ++k)
    for (Int i = 0; i < len; ++i)
      (*density_p[k])[i] = 1;
  std::vector<Real> wrk(opts.d2c ? 2*dnn : cnn);
  // Record the initial and final states.
  std::vector<Real> error_data(4*len);

  {
    // Get the initial conditions.
    std::vector<Real> lat(cnn), lon(cnn);
    for (Int i = 0; i < cnn; ++i) {
      const auto n = slice(m.cgll_p, i);
      xyz2ll(n[0], n[1], n[2], lat[i], lon[i]);
    }
    Real* data = opts.d2c ? tracer_p[0]->data() : wrk.data();
    gallery::InitialCondition::init(
      ic, nslices(m.cgll_p), lat.data(), lon.data(), data);
    // Record the ICs.
    if ( ! opts.d2c)
      map_cgll2dgll(m.dglln2cglln, data, tracer_p[0]->data());
    memcpy(error_data.data(), tracer_p[0]->data(), len*sizeof(Real));
    if (ncw) {
      ncw->advance_time_to(0);
      ncw->write_field("tracer", data);
    }
    memcpy(error_data.data() + len, density_p[0]->data(), len*sizeof(Real));
    if (ncw) {
      data = opts.d2c ? density_p[0]->data() : wrk.data();
      if ( ! opts.d2c)
        map_dgll2cgll(m.dglln2cglln, rd.dgbfi(), rd.cgbfi(),
                      density_p[0]->data(), data, cnn);
      ncw->write_field("density", data);
    }
  }
  // Remap is done on density*tracer, but sometimes the tracer field doesn't
  // have the density rho in it.
  for (Int i = 0; i < len; ++i)
    (*tracer_p[0])[i] *= (*density_p[0])[i];

  // Time step.
  Vec3s::HostMirror departure_p;
  ko::resize(departure_p, nslices(m.geo_p), szslice(m.geo_p));
  const Real dt = T/nsteps;
  const Int last_step =
    opts.stepping == IntegrateOptions::test_looa ? 1 : nsteps - 1;
  ProgressBar progress_bar("integrate", last_step+1, 10);
  const auto step_t = siqk::tic();
  Timer::stop(Timer::ts_setup); Timer::start(Timer::ts);
  for (Int step = 0; step <= last_step; ++step) {
    Timer::start(Timer::ts_integrate);
    const Real tf = step == last_step ? T : dt*(step + 1);
    switch (opts.stepping) {
    case IntegrateOptions::fwd:
      // Integrate mesh forward in time.
      mi->integrate(dt*step, tf, departure_p);
      break;
    case IntegrateOptions::bwd:
      throw std::runtime_error("IntegrateOptions::bwd is not impl'ed.");
      break;
    case IntegrateOptions::test_looa:
      switch (step) {
      case 0: mi->integrate(dt*step, tf, departure_p); break;
      case 1: mi->integrate(dt*step, dt*(step - 1), departure_p); break;
      default: assert(0); break;
      }
      break;
    }
    Timer::stop(Timer::ts_integrate); Timer::start(Timer::ts_remap);
    remap(rd, m, departure_p, tracer_p[0]->data(), tracer_p[1]->data(), 1,
          density_p[0]->data(), density_p[1]->data(), ! opts.d2c, wrk.data());
    Timer::stop(Timer::ts_remap); Timer::start(Timer::ts_rest);
    if (step == 0) {
      // Analyze the remap operator R = M \ T.
      RealArray::HostMirror dgbfi_s;
      calc_basis_function_integrals(m.np, m.tq_order, departure_p, m.geo_c2n,
                                    dgbfi_s);
      printf("\n> triangle quadrature jacobians\n");
      rd.check(dgbfi_s.ptr_on_device(), rd.dgbfi().ptr_on_device());
      RealArray::HostMirror Js_gll;
      calc_gll_basis_function_integrals(m, departure_p, Js_gll);
      printf("> GLL jacobians\n");
      rd.check(Js_gll.ptr_on_device(), Jt_gll.ptr_on_device());
    }
    gdbg.write("T", rd.T());
    gdbg.write_p("geo_p", m.geo_p); gdbg.write_c2n("geo_c2n", m.geo_c2n);
    gdbg.write_p("departure_p", departure_p);

    // Netcdf I/O.
    if (ncw && (step % write_every == 0 || step == last_step)) {
      ncw->advance_time_to(tf);
      if (opts.d2c) {
        ncw->write_field("tracer", tracer_p[1]->data());
        ncw->write_field("density", density_p[1]->data());
      } else {
        Real* const data = wrk.data();
        map_dgll2cgll(m.dglln2cglln, rd.dgbfi(), rd.cgbfi(),
                      tracer_p[1]->data(), data, cnn);
        ncw->write_field("tracer", data);
        map_dgll2cgll(m.dglln2cglln, rd.dgbfi(), rd.cgbfi(),
                      density_p[1]->data(), data, cnn);
        ncw->write_field("density", data);
      }
    }
    // Record data for error analysis.
    if (step == last_step) {
      memcpy(error_data.data() + 2*len, tracer_p[1]->data(), len*sizeof(Real));
      memcpy(error_data.data() + 3*len, density_p[1]->data(), len*sizeof(Real));
    }

#   pragma omp parallel for
    for (Int i = 0; i < len; ++i)
      (*tracer_p[1])[i] *= (*density_p[1])[i];

    std::swap(tracer_p[0], tracer_p[1]);
    std::swap(density_p[0], density_p[1]);
    progress_bar.update();
    gdbg.advance();
    gdbg.set_on(false);
  }
  const Real step_et = siqk::toc(step_t);
  Timer::stop(Timer::ts);
  siqk::print_times("timestep", step_et);
  out.et_timestep = step_et;

  { Timer::start(Timer::ts_error);
    const Real* const d = error_data.data();
    print_error(m, Jt_gll, ! opts.d2c, d, d + len, d + 2*len,
                d + 3*len, out);
    Timer::stop(Timer::ts_error);
  }
}

static void run (const Input& in) {
  const Real T = day2sec(in.T);
  const auto ic = gallery::InitialCondition::from_string(
    in.initial_condition);
  RemapOptions ro;
  ro.np = in.np;

  Mesh m;
  init_mesh(in.np, in.tq_order, in.ne, m);

  auto mi = MeshIntegratorFactory::create(in.ode, in.xyz_form, m.geo_p);
  // Get lat-lon of geo mesh nodes.
  const Int nn = nslices(m.geo_p);
# pragma omp parallel for
  for (Int i = 0; i < nn; ++i) {
    const auto n = slice(m.geo_p, i);
    Real* const lli = mi->get_ll().data() + 2*i;
    xyz2ll(n[0], n[1], n[2], lli[0], lli[1]);
  }

  Output out;
  integrate(m, mi, ro, T, in.nsteps, ic, in.output_fn, in.write_every,
            in.integrate_options, out);
  print_one_liner(in, out);
}

Input::Input (int argc, char** argv)
  : output_fn("tmp/out"), ode("divergent"),
    initial_condition("xyztrig"), T(12), ne(5), nsteps(120), write_every(1),
    monotone_type(0), np(4), tq_order(12), debug(false), xyz_form(false)
{
  program_name = argv[0];
  integrate_options.stepping = IntegrateOptions::fwd;
  integrate_options.d2c = false;
  for (int i = 1; i < argc; ++i) {
    const std::string& token = argv[i];
    if (eq(token, "-o", "--output"))
      output_fn = argv[++i];
    else if (eq(token, "-T"))
      T = atof(argv[++i]);
    else if (eq(token, "-nsteps"))
      nsteps = atoi(argv[++i]);
    else if (eq(token, "-ode"))
      ode = argv[++i];
    else if (eq(token, "-ic"))
      initial_condition = argv[++i];
    else if (eq(token, "-mono", "--monotone"))
      monotone_type = atoi(argv[++i]);
    else if (eq(token, "-np"))
      np = atoi(argv[++i]);
    else if (eq(token, "-tq"))
      tq_order = atoi(argv[++i]);
    else if (eq(token, "-ne"))
      ne = atoi(argv[++i]);
    else if (eq(token, "-we", "--write-every"))
      write_every = atoi(argv[++i]);
    else if (eq(token, "-looa", "--looa"))
      integrate_options.stepping = IntegrateOptions::test_looa;
    else if (eq(token, "-xyz", "--xyz"))
      xyz_form = true;
    else if (eq(token, "-d2c", "--d2c"))
      integrate_options.d2c = true;
    else if (eq(token, "-d", "--debug"))
      debug = true;
  }

  if (np == 4) tq_order = 20;

  print(std::cout);
}

void Input::print (std::ostream& os) const {
  os << "output filename (-o): " << output_fn << "\n"
     << "ode (-ode, " << MeshIntegratorFactory::get_inputs() << "): "
     <<   ode << "\n"
     << "xyz_form (-xyz): " << xyz_form << "\n"
     << "initial condition (-ic, "
     <<   gallery::InitialCondition::get_inputs() << "): "
     <<   initial_condition << "\n"
     << "T (-T): " << T << " [day]\n"
     << "nsteps (-nsteps): " << nsteps << "\n"
     << "np (-np): " << np << "\n"
     << "tq (-tq): " << tq_order << "\n"
     << "ne (-ne): " << ne << "\n"
     << "monotone_type (-mono, {0,1,2,3}): " << monotone_type << "\n"
     << "write every (-we): " << write_every << "\n"
     << "test_looa (--looa): "
     << (integrate_options.stepping == IntegrateOptions::test_looa) << "\n"
     << "d2c (-d2c): " << integrate_options.d2c << "\n"
     << "debug (-d): " << debug << "\n";
}

int main (int argc, char** argv) {
  Kokkos::initialize(argc, argv); {
    Timer::init();
    Timer::start(Timer::total);
    BlockMatrix<>::test();
    Input in(argc, argv);
    run(in);
    Timer::stop(Timer::total);
    Timer::print();
  } Kokkos::finalize_all();
}
