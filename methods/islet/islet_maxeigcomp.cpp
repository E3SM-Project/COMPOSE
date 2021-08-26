#include <omp.h>

#include "islet_maxeigcomp.hpp"
#include "islet_tables.hpp"
#include "islet_npx.hpp"
#include "islet_util.hpp"

// LAPACK eigendecomp routine for real unsymmetric matrix.
typedef int fint;
extern "C" void dgeev_ (char* jobvl, char* jobvr, fint* n, double* a, int* lda,
                        double* wr, double* wi,
                        double* vl, int* ldvl,
                        double* vr, int* ldvr,
                        double* work, int* lwork, int* info);
extern "C" void zgeev_ (char* jobvl, char* jobvr, fint* n, Complex* a, int* lda,
                        Complex* w,
                        Complex* vl, int* ldvl,
                        Complex* vr, int* ldvr,
                        Complex* work, int* lwork, double* rwork, int* info);
// LAPACK SVD routine for real unsymmetric matrix.
extern "C" void dgesvd_ (char* jobu, char* jobvt, fint* m, fint* n, double* a, int* lda,
                         double* s, double* u, int* ldu, double* vt, int* ldvt,
                         double* work, int* lwork, int* info);
extern "C" void zgesvd_ (char* jobu, char* jobvt, fint* m, fint* n, Complex* a, int* lda,
                         double* s, Complex* u, int* ldu, Complex* vt, int* ldvt,
                         Complex* work, int* lwork, double* rwork, int* info);

static
void dgeev (char jobvl, char jobvr, int n, double* a, int lda,
            double* wr, double* wi,
            double* vl, int ldvl,
            double* vr, int ldvr,
            double* work, int lwork, int& info) {
  dgeev_(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
         work, &lwork, &info);
}

static
void zgeev (char jobvl, char jobvr, int n, Complex* a, int lda,
            Complex* w,
            Complex* vl, int ldvl,
            Complex* vr, int ldvr,
            Complex* work, int lwork, double* rwork, int& info) {
  zgeev_(&jobvl, &jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr,
         work, &lwork, rwork, &info);
}

static
void dgesvd (char jobu, char jobvt, fint m, fint n, double* a, int lda,
             double* s, double* u, int ldu, double* vt, int ldvt,
             double* work, int lwork, int& info) {
  dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
          work, &lwork, &info);
}

static
void zgesvd (char jobu, char jobvt, fint m, fint n, Complex* a, int lda,
             double* s, Complex* u, int ldu, Complex* vt, int ldvt,
             Complex* work, int lwork, double* rwork, int& info) {
  zgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
          work, &lwork, rwork, &info);
}

class Mesh {
  Int nc_;
  Real dx_;
public:
  Mesh (const Int nc) { init(nc); }
  void init (const Int nc) {
    nc_ = nc;
    dx_ = 1.0/nc;
  }

  Real toperiodic (const Real& x) const { return x - std::floor(x); }
  Int ncell () const { return nc_; }
  Real dx () const { return dx_; }
  Int incell (const Real& x) const { return std::floor(toperiodic(x) / dx_); }
  void toref (const Real& x, Int& ci, Real& a) const {
    ci = incell(x);
    a = 2*(toperiodic(x)*nc_ - ci) - 1;
  }
  Real tophysical (const Int& ci, const Real& a) const {
    return toperiodic((ci + 0.5*(a+1))*dx_);
  }

  static Int unittest () {
    Mesh m(42);
    const Real dx = m.dx();
    const Int nc = m.ncell();
    const Real eps = std::numeric_limits<Real>::epsilon();
    using islet::reldif;
    Int ne = 0;
    if (reldif(m.toperiodic( 1.7), 0.7) > 10*eps) ++ne;
    if (reldif(m.toperiodic(-0.8), 0.2) > 10*eps) ++ne;
    if (m.incell(2) != 0) ++ne;
    if (m.incell(2 - 0.5*dx) != nc-1) ++ne;
    { Int ic; Real a;
      const Real x = 4.4*dx;
      m.toref(x, ic, a);
      if (ic != 4) ++ne;
      if (reldif(a, -0.2) > 1e2*eps) ++ne;
      if (reldif(m.tophysical(ic, a), x) > 1e2*eps) ++ne; }
    return ne;
  }
};

void op_apply (const Mesh& m, const Int ne, const InterpMethod& im,
               const Real dx_flow, const Real* const src, Real* const tgt) {
  const auto xnodes = im.get_xnodes();
  for (Int ci = 0, k = 0; ci < ne; ++ci) {
    for (Int i = 0; i < im.np-1; ++i, ++k) {
      Int ci_src;
      Real a_src;
      m.toref(m.tophysical(ci, xnodes[i]) + dx_flow, ci_src, a_src);
      Real v[32];
      op_eval(im, a_src, v);
      Real val = 0;
      const Real* src_cell = src + (im.np - 1)*ci_src;
      for (Int i_src = 0; i_src < im.np; ++i_src)
        val += v[i_src]*src_cell[i_src];
      tgt[k] = val;
    }
  }
}

static void get_matrix (
  const Int& ne, const Int& np, const Real& dx, const InterpMethod& interp_method,
  Array<Real>& A)
{
  // Get the CSL operator as a matrix.
  const Int N = ne*(np-1), Np1 = N+1;
  A.optclear_and_resize(N*N);
  Mesh m(ne);
  const Real dx_flow = -dx/ne;
  Array<Real> u0(Np1, 0);
  const Real* const src = u0.data();
  for (Int b = 0; b < N; ++b) {
    // Get b'th column.
    u0[b] = 1;
    if (b == 0) u0[N] = 1;
    else if (b == N) u0[0] = 1;
    Real* const tgt = A.data() + b*N;
    op_apply(m, ne, interp_method, dx_flow, src, tgt);
    u0[b] = 0;
    if (b == 0) u0[N] = 0;
    else if (b == N) u0[0] = 0;
  }
}   

// A is overwritten.
static Real cond_2norm (Complex* A, const Int& n, Real* work, const Int& nwork) {
  throw_if(nwork < 12*n, "work should be >= 12 n");
  Real* rwork = work;
  Real* s = work + 5*n;
  Complex* cwork = reinterpret_cast<Complex*>(work + 6*n);
  const Int lwork = (nwork - 6*n)/2;
  int info;
  zgesvd('n', 'n', n, n, A, n, s,
         nullptr, 1, nullptr, 1,
         cwork, lwork, rwork,
         info);
  return s[0]/s[n-1];
}

namespace bloch {
// K is the number of nodes still in the element after shifting right by the
// fraction dx of an element.
static Int get_K (const Int np, const Real* x_gll, const Real dx) {
  assert(dx > 0 && dx < 1);
  Int K;
  for (K = 0; K < np; ++K)
    if (x_gll[K] + 2*dx >= 1)
      break;
  assert(x_gll[0] > -1 || (K > 0 && K < np));
  return K;
}

// K is the opposite of what it is in csl.hy.
struct Data {
  const Int np;
  const Real* const x_gll;
  const Int K;
  const Real dx;

  Data (const Int& inp, const Real& idx, const Real* const x_gll_ = nullptr)
    : np(inp),
      x_gll(x_gll_ ? x_gll_ : islet::get_x_gll(np)),
      K(get_K(inp, x_gll, idx)),
      dx(idx)
  {}
};

// Get the (np-1)xnp block matrix that is repeated in A.
void form_kernel_block (const Data& d, const InterpMethod& im, Real* const A) {
  const Int npm1 = d.np-1;
  // Get the kernel block in an (np-1)xnp row-major matrix starting at A.
  for (Int ip = 0; ip < npm1; ++ip) {
    Real ref = d.x_gll[ip] + 2*d.dx;
    if (ref >= 1) {
      assert(ip >= d.K);
      ref -= 2;
    } else {
      if (ip >= d.K) {
        pr(puf(d.np) pu(d.K) pu(d.dx) pu(ip) pu(ref));
        islet::prarr("d.x_gll", d.x_gll, d.np);
      }
      assert(ip < d.K);
    }
    assert(ref >= -1 && ref < 1);
    op_eval(im, ref, A + d.np*ip);
  }
}

// Get B(mu) from A.
void form_Bmu (const Data& d, const Complex& mu, const Real* const A,
               Complex* const Bmu) {
  const Int npm1 = d.np-1;
  for (Int c = 0; c < npm1; ++c) {
    Complex* const col = Bmu + npm1*c;
    for (Int r = 0  ; r < d.K ; ++r) col[r] =    A[d.np*r + c];
    for (Int r = d.K; r < npm1; ++r) col[r] = mu*A[d.np*r + c];
  }
  {
    Complex* const col = Bmu;
    for (Int r = 0  ; r < d.K ; ++r) col[r] +=    mu*A[d.np*r + npm1];
    for (Int r = d.K; r < npm1; ++r) col[r] += mu*mu*A[d.np*r + npm1];
  }  
}

void form_Bmu (const Data& d, const Int& ne, const Int& ie, const Real* const A,
               Complex* const Bmu) {
  const Real arg = 2 * M_PI * (Real(ie)/ne);
  const Complex mu(std::cos(arg), std::sin(arg));
  form_Bmu(d, mu, A, Bmu);
}

void edecomp (const Data& d, Complex* const Bmu,
              Real* const work, Int lwork,
              // (w[2*i], ws(2*i+1)) contains the i'th eigenvalue.
              Real* const w, Complex* const V = nullptr) {
  const Int npm1 = d.np-1;
  assert(lwork >= 22*npm1);
  lwork -= 2*npm1;
  Int eig_info;
  zgeev('n', V ? 'v' : 'n', npm1, Bmu, npm1,
        reinterpret_cast<Complex* const>(w),
        nullptr, 1, V, npm1,
        reinterpret_cast<Complex* const>(work), lwork/2,
        work + lwork, eig_info);
  assert(eig_info == 0);
}
} // namespace bloch

void MaxEigComputer::setup_workspace (const Int max_ne_) {
  max_ne = max_ne_;
  const int nthr = threaded ? omp_get_max_threads() : 1;
  wss.resize(nthr);
}

struct IndexCover {
  IndexCover (Int n_, Int P_) {
    n = n_;
    P = P_;
    split = std::max(10, P);
    N = std::max(1, n/split);
  }

  Int nit () const { return (n+P-1)/P; }

  Int idx (const Int it, const Int tid) const {
    const Int k = P*it + tid;
    const Int i = (k >= split*N) ? k : (k % split)*N + (k / split);
    return i >= n ? -1 : i;
  }

private:
  Int P, n, N, split;
};

static Int test_index_cover () {
  Int nerr = 0;
  const auto check = [&] (const Int n, const Int P) -> Int {
    Int ne = 0;
    std::vector<int> cnt(n, 0);
    IndexCover ic(n, P);
    const Int nit = ic.nit();
    for (Int it = 0; it < nit; ++it)
      for (Int tid = 0; tid < P; ++tid) {
        const Int i = ic.idx(it, tid);
        if (i >= n) ++ne;
        else if (i >= 0) ++cnt[i];
      }
    for (Int i = 0; i < n; ++i)
      if (cnt[i] != 1)
        ++ne;
    return ne;
  };
  for (const  Int n : {15, 33, 128, 1111, 3333, 4000, 7777})
    for (const Int P : {1, 2, 3, 8, 11, 48, 272})
      nerr += check(n, P);
  return nerr;
}

Real MaxEigComputer::
run (const Int& ne_max, const Int& ndx_max, const Real& maxeigampm1,
     const bool quiet, const InterpMethod& im) {
  setup_workspace(ne_max);
  // Search dx in (0, 0.5], in parallel, to see if there's a max |lambda| - 1
  // bigger than tol.
  const auto cdxeig = [=] (Int ne, Int ndx, Real tol) {
    Real mme = 0;
    const int P = threaded ? omp_get_max_threads() : 1;
    // Chunk up the search space so we explore widely as early as possible.
    const bool both_dir = false;
    const auto fac = both_dir ? 2 : 1;
    const Int n = fac*ndx;
    IndexCover ic(n, P);
    const Int nit = ic.nit();
    assert(P <= max_nthread);
    const auto run1 = [&] (const int it, const int tid) -> Real {
      const auto i = ic.idx(it, tid) + 1;
      if (i < 0) return 1;
      // dx is in (-0.5, 0.5] or (0, 0.5].
      const Real dx = both_dir ? (Real(i)/n - 0.5) : (0.5*i)/n;
      if (dx == 0) return 1;
      Real me;
      compute(im, dx, ne, &me);
      return me;
    };
    std::array<Real,max_nthread> mes;
    for (int it = 0; it < nit; ++it) {
      if (threaded) {
#       pragma omp parallel
        {
          const int tid = omp_get_thread_num();
          mes[tid] = run1(it, tid);
        }
        for (int j = 0; j < std::min(P, n); ++j)
          mme = std::max(mme, mes[j]);
      } else {
        mme = run1(it, 0);
      }
      if (mme - 1 >= tol)
        break;
    }
    return mme;
  };

  Real maxeigamp = -1;
  // Ramp up precision of search to encourage an early exit when the method
  // (np, order, offset) is not stable.
  for (int ne : {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, ne_max}) {
    if (ne > ne_max) break;
    bool toobig = false, first = true;
    for (int ndx : {32, 64, 128, 256, 512, 1024, ndx_max}) {
      if ( ! first && ndx > ndx_max) break;
      if (ne < ne_max && ndx == ndx_max) break;
      first = false;
      maxeigamp = std::max(maxeigamp,
                           cdxeig(ne, ndx, maxeigampm1));
      if ( ! quiet)
        printf("ne %2d nsamp %3d %1.3e\n", ne, ndx, maxeigamp-1);
      toobig = maxeigamp - 1 > maxeigampm1;
      if (toobig) break;
    }
    if (toobig) break;
  }
  return maxeigamp-1;
}

MaxEigComputer::Analysis MaxEigComputer::
calc_max_vals (const Int& nmu, const Int& ndx,
               const InterpMethod& im) {
  const auto calc1 = [=] (const Int& idx) -> Analysis {
    const Real dx = (0.5*idx)/ndx;
    if (dx == 0) return {1,1,1};
    Analysis vals;
    compute(im, dx, nmu, &vals.max_eig_amp, &vals.max_condv,
            &vals.max_defect_ub);
    return vals;
  };
  Real max_eig_amp = 0, max_condv = 0;
  // Use <= ndx below to get dx in [0, 0.5] instead of [0, 0.5).
  if (threaded) {
    Real mea[max_nthread] = {0}, mcv[max_nthread] = {0};
#   pragma omp parallel for
    for (Int ix = 0; ix <= ndx; ++ix) {
      const int tid = omp_get_thread_num();
      const auto mv = calc1(ix);
      mea[tid] = std::max(mea[tid], mv.max_eig_amp);
      mcv[tid] = std::max(mcv[tid], mv.max_condv);
    }
    for (Int i = 0; i < omp_get_max_threads(); ++i) {
      max_eig_amp = std::max(max_eig_amp, mea[i]);
      max_condv = std::max(max_condv, mcv[i]);
    }
  } else {
    for (Int ix = 0; ix <= ndx; ++ix) {
      const auto mv = calc1(ix);
      max_eig_amp = std::max(max_eig_amp, mv.max_eig_amp);
      max_condv = std::max(max_condv, mv.max_condv);
    }
  }
  return {max_eig_amp, max_condv, 1};
}

void MaxEigComputer::
compute (const InterpMethod& im, const Real& dx, const Int& ne,
         Real* max_amp_out, Real* max_condv,
         Real* max_defect_ub,
         Complex* lam, Complex* V) {
  auto& ws = wss[threaded ? omp_get_thread_num() : 0];
  Real max_amp = 0;
  const Int np = im.np, npm1 = np - 1;
  if (bloch) {
    const Int N = 3*npm1*np, edecomp_ws = 22*npm1, evecs_ws = 2*npm1*npm1,
      cond_2norm_ws = 12*npm1;
    if (static_cast<int>(ws.A.size()) < N) {
      ws.A.optclear_and_resize(N);
      ws.wr.optclear_and_resize(2*npm1);
      ws.work.optclear_and_resize(edecomp_ws + evecs_ws + cond_2norm_ws);
    }
    bloch::Data bd(np, dx, im.uim == nullptr ? nullptr : im.uim->get_xnodes());
    bloch::form_kernel_block(bd, im, ws.A.data());
    Complex* const Bmu = reinterpret_cast<Complex*>(ws.A.data() + npm1*np);
    if (max_condv) *max_condv = 0;
    Complex* evecs = reinterpret_cast<Complex*>(ws.work.data() + edecomp_ws);
    for (int ie = 0; ie < ne; ++ie) {
      bloch::form_Bmu(bd, ne, ie, ws.A.data(), Bmu);
      bloch::edecomp(bd, Bmu, ws.work.data(), edecomp_ws, ws.wr.data(),
                     V || max_condv ? evecs : nullptr);
      for (Int i = 0; i < npm1; ++i) {
        const Real re = ws.wr[2*i], im = ws.wr[2*i+1];
        max_amp = std::max(max_amp, std::sqrt(re*re + im*im));
      }
      if (V)
        for (Int i = 0; i < npm1*npm1; ++i)
          V[npm1*npm1*ie + i] = evecs[i];
      if (max_condv) {
        const Real condv = cond_2norm(evecs, npm1,
                                      ws.work.data() + edecomp_ws + evecs_ws,
                                      cond_2norm_ws);
        *max_condv = std::max(*max_condv, condv);
      }
      if (lam)
        for (Int i = 0; i < npm1; ++i)
          lam[ie*npm1 + i] = Complex(ws.wr[2*i], ws.wr[2*i+1]);
    }
  } else {
    throw_if(max_condv || max_defect_ub,
             "Only Bloch-wave-based edecomp supports"
             " max cond(V) and max defect_ub(lam).");
    // Compute eigenvalues of discrete space-time operator with ne cells.
    const Int N = ne*npm1;
    if (static_cast<int>(ws.wr.size()) < N) {
      const int Nonce = std::max(ne, max_ne)*npm1;
      ws.A.optclear_and_resize(Nonce*Nonce);
      ws.wr.optclear_and_resize(Nonce);
      ws.wi.optclear_and_resize(Nonce);
      ws.work.optclear_and_resize(10*Nonce);
    }
    // -dx to match Bloch-wave computation.
    get_matrix(ne, np, -dx, im, ws.A);
    int eig_info;
    Real* Vreal = V ? reinterpret_cast<Real*>(V) : nullptr;
    dgeev('n', V ? 'v' : 'n', N, ws.A.data(), N, ws.wr.data(), ws.wi.data(),
          nullptr, 1, Vreal, N, ws.work.data(), ws.work.size(), eig_info);
    for (Int i = 0; i < N; ++i)
      max_amp = std::max(max_amp,
                         std::sqrt(ws.wr[i]*ws.wr[i] + ws.wi[i]*ws.wi[i]));
    if (lam)
      for (Int i = 0; i < N; ++i)
        lam[i] = Complex(ws.wr[i], ws.wi[i]);
  }
  assert(max_amp_out);
  *max_amp_out = max_amp;
}

static void normalize (Complex* v, const Int& n) {
  Real norm2 = 0;
  for (Int k = 0; k < n; ++k) {
    Real d = std::abs(v[k]);
    norm2 += d*d;
  }
  const Real scale = 1/std::sqrt(norm2);
  for (Int k = 0; k < n; ++k)
    v[k] *= scale;
}

static void remove_phase (Complex* v, const Int& n) {
  const Real arg = std::atan2(v[0].imag(), v[0].real());
  Complex phase(std::cos(arg), -std::sin(arg));
  for (Int k = 0; k < n; ++k)
    v[k] *= phase;
}

static void write_dgeev_evec (const Complex& lam, const Real* V,
                              const Int& vi, const Int& pair, const Int& n,
                              Complex* v) {
  if (lam.imag() == 0) {
    for (Int k = 0; k < n; ++k)
      v[k] = V[vi*n + k];
  } else {
    const Int sign = pair == 0 ? 1 : -1;
    for (Int k = 0; k < n; ++k)
      v[k] = Complex(V[vi*n + k], sign*V[(vi+1)*n + k]);
  }
  normalize(v, n);
  remove_phase(v, n);
}

static void write_bloch_evec (const Complex* v_bloch, const Int& vi,
                              const Int& ne, const Int& np, Complex* v) {
  const Int npm1 = np-1, n = ne*npm1;
  const Int vie = vi / npm1;
  const Real arg = 2 * M_PI * (Real(vie)/ne);
  const Complex* const phi = v_bloch + npm1*vi;
  for (Int ie = 0; ie < ne; ++ie) {
    const Complex mu(std::cos(ie*arg), std::sin(ie*arg));
    for (Int i = 0; i < npm1; ++i)
      v[ie*npm1 + i] = mu*phi[i];
  }
  normalize(v, n);
  remove_phase(v, n);    
}

// Check that the eigenvalues and eigenvectors derived from Bloch-wave analysis
// -- the size-(np-1) eigenvalue problem -- match those derived from brute-force
// eigendecomp of the ne-mesh space-time operator.
static int
check_bloch_edecomp (const int np, const int ne,
                     const Complex* lam, const Complex* lam_bloch,
                     const Real* V, const Complex* v_bloch,
                     // >= 2 ne (np-1)
                     Complex* work) {
  const Int N = ne*(np-1);
  Complex* const v_bf = work;
  Complex* const v_full_bloch = work + N;
  int nerr = 0;
  Array<bool> used(N, false);
  for (Int i = 0, i_vec = 0, pair = 0; i < N; ++i) {
    // Find corresponding eigenvalues.
    const auto lam_brute_force = lam[i];
    Real min_diff = 2;
    Int jmin = -1;
    for (Int j = 0; j < N; ++j) {
      const Real diff = std::abs(lam_bloch[j] - lam_brute_force);
      if (diff < min_diff) {
        min_diff = diff;
        jmin = j;
      }
    }
    used[jmin] = true;
    if (min_diff > 1e-13)
      ++nerr;
    // Compare eigenvectors. For dgeev, some arithmetic is needed to account
    // for how evecs are packed.
    write_dgeev_evec(lam[i], V, i_vec, pair, N, v_bf);
    if (pair == 1) {
      assert(lam[i].imag() != 0);
      i_vec += 2;
      pair = 0;
    } else if (lam[i].imag() == 0) {
      i_vec++;
      assert(pair == 0);
    } else {
      assert(pair == 0);
      pair++;
    }
    write_bloch_evec(v_bloch, jmin, ne, np, v_full_bloch);
    Real num = 0, den = 0;
    for (Int k = 0; k < N; ++k) {
      Real d = std::abs(v_full_bloch[k] - v_bf[k]);
      num += d*d;
      d = std::abs(v_bf[k]);
      den += d*d;
    }
    if (std::sqrt(num/den) >= 3e-9) {
      pr(puf(i) pu(jmin) pu(min_diff) pu(std::sqrt(num/den)));
      ++nerr;
    }
  }
  // Check that all eigenvalues were matched.
  for (Int j = 0; j < N; ++j)
    if ( ! used[j])
      ++nerr;
  return nerr;
}

int MaxEigComputer::unittest () {
  int nerr = 0;
  MaxEigComputer mec(false, false), mec_bloch(false, true);
  const int np_max = std::min(13, islet::np_max), ne_max = np_max,
    Nmax = (np_max-1)*ne_max;
  Array<Complex> lam(Nmax), lam_bloch(Nmax),
    V((Nmax*Nmax)/2), v_bloch((np_max-1)*(np_max-1)*ne_max), work(2*Nmax);
  for (int np = 4; np <= np_max; ++np) {
    InterpMethod im(np, InterpMethod::npxstab);
    for (int ne : {1, 3, ne_max}) {
      const Int N = ne*(np-1);
      assert(N <= Nmax);
      for (Real dx : {0.05, 0.42, 0.7}) {
        Real max_amp;
        mec.compute(im, dx, ne, &max_amp, nullptr, nullptr,
                    lam.data(), V.data());
        mec_bloch.compute(im, dx, ne, &max_amp, nullptr, nullptr,
                          lam_bloch.data(), v_bloch.data());
        nerr += check_bloch_edecomp(np, ne, lam.data(), lam_bloch.data(),
                                    reinterpret_cast<Real*>(V.data()),
                                    v_bloch.data(), work.data());
      }
    }
  }
  nerr += test_index_cover();
  return nerr;
}
