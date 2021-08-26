/* This standalone program implements the example unstable 1D and 2D problems
   that use the classical cubic interpolation semi-Lagrangian method.

   Build it:
     g++ -O3 -c cslunstab.cpp
     g++ cslunstab.o -llapack -lblas -o cslunstab
   Run it:
     ./cslunstab
   There will be no output if all the assertions pass. The two lines
     require(cubic1d_demo_unstable_problem() >= 1 + 1e-3);
   and
     require(cubic2d_demo_unstable_problem() >= 1 + 1e-2);
   assert that the maximum eigenvalue magnitude is at least 1 + a small amount
   in each problem.
 */

#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <sstream>

using Int = int;
using Real = double;

#define require(condition) do {                                         \
    if ( ! (condition)) {                                               \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": FAIL:\n" << #condition  \
           << "\n";                                                     \
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

bool eq(const std::string& a, const char* const b1, const char* const b2 = 0);

template <typename T> T square (const T& x) { return x*x; }

inline Real reldif (const Real& a, const Real& b)
{ return std::abs(b - a)/std::max(std::abs(a), std::abs(b)); }

template <typename T1, typename T2>
inline bool equal (const T1& a, const T2& b) {
  if (a != b)
    printf("equal: a,b = %23.16e %23.16e re = %23.16e\n",
           Real(a), Real(b), std::abs((a-b)/Real(a)));
  return a == b;
}

template <typename T1, typename T2>
inline bool almost_equal (const T1& a, const T2& b, const Real tol) {
  const auto re = std::abs(a-b)/(1.0 + std::abs(a));
  const bool good = re <= tol;
  if ( ! good)
    printf("equal: a,b = %23.16e %23.16e re = %23.16e tol %9.2e\n",
           Real(a), Real(b), re, tol);
  return good;
}

inline double urand () { return rand() / ((double) RAND_MAX + 1.0); }

extern "C" void dgeev_(char* jobvl, char* jobvr, int* n, double* a, int* lda,
                       double* wr, double* wi,
                       double* vl, int* ldvl,
                       double* vr, int* ldvr,
                       double* work, int* lwork, int* info);

void dgeev (int n, double* a, int lda,
            double* wr, double* wi,
            std::vector<double>& work, int& info,
            double* vl = nullptr, int ldvl = 1,
            double* vr = nullptr, int ldvr = 1) {
  int lwork = 10*n;
  if (static_cast<int>(work.size()) < lwork) work.resize(lwork);
  char jobvl = vl ? 'v' : 'n';
  char jobvr = vr ? 'v' : 'n';
  assert(vl == nullptr || (ldvl >= n));
  assert(vr == nullptr || (ldvr >= n));
  dgeev_(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
         work.data(), &lwork, &info);
}

Real calc_max_eig_amp (int n, double* a, int lda, std::vector<double>& work) {
  work.resize(12*n);
  Real* const wr = work.data() + 10*n;
  Real* const wi = work.data() + 11*n;
  int info;
  dgeev(n, a, lda, wr, wi, work, info);
  Real mea = 0;
  for (int i = 0; i < n; ++i) mea = std::max(mea, std::sqrt(square(wr[i]) + square(wi[i])));
  return mea;
}

static void eval_lagrange_poly_basis (const Int& n, const Real* xsup, const Real& x,
                                      Real* const y) {
  for (int i = 0; i < n; ++i) {
    Real f = 1;
    for (int j = 0; j < n; ++j)
      f *= (i == j) ?
        1 :
        (x - xsup[j]) / (xsup[i] - xsup[j]);
    y[i] = f;
  }
}

static Real eval_lagrange_poly (const Int& n, const Real* xsup, const Real* ysup,
                                const Real& x) {
  Real y = 0;
  for (int i = 0; i < n; ++i) {
    Real f = 1;
    for (int j = 0; j < n; ++j)
      f *= (i == j) ?
        1 :
        (x - xsup[j]) / (xsup[i] - xsup[j]);
    y += f*ysup[i];
  }
  return y;
}

// Move x to its in-bounds periodic point.
static Real calc_periodic_x (const Real& x0, const Real& x1, const Real x) {
  if (x >= x0 && x <= x1) return x;
  const auto D = x1 - x0;
  auto xper = (x - x0)/D;
  xper = x0 + (xper - std::floor(xper))*D;
  assert(xper >= x0);
  assert(xper <= x1);
  return xper;
}

// Find xper's cell.
static Int find_cell (const Real* const xnode, const Int nreg, const Real& xper) {
  if (xper == xnode[0]) return 0;
  const Int iper = static_cast<Int>(
    std::lower_bound(xnode, xnode + nreg + 1, xper)
    - xnode)
    - 1;
  assert(iper >= 0   );
  assert(iper <  nreg);
  assert(xper >= xnode[iper  ]);
  assert(xper <= xnode[iper+1]);
  return iper;
}

// xnode[0] is identical to xnode[nreg]. ynode[0] must equal ynode[nreg]. We
// take ynode[nreg] as input rather than use ynode[0] to be sure the caller has
// set up a periodic problem.
static Real periodic_cubic1d_interp (
  const Int nreg, const Real* const xnode, const Real* const ynode, const Real& xeval)
{
  assert(ynode[0] == ynode[nreg]);

  const auto xper = calc_periodic_x(xnode[0], xnode[nreg], xeval);
  const auto iper = find_cell(xnode, nreg, xper);

  Real xsup[4], ysup[4];
  if (iper == 0) {
    xsup[0] = xnode[0] - (xnode[nreg] - xnode[nreg-1]);
    ysup[0] = ynode[nreg-1];
    for (int i = 1; i < 4; ++i) xsup[i] = xnode[iper-1+i];
    for (int i = 1; i < 4; ++i) ysup[i] = ynode[iper-1+i];
  } else if (iper == nreg-1) {
    for (int i = 0; i < 3; ++i) xsup[i] = xnode[iper-1+i];
    for (int i = 0; i < 3; ++i) ysup[i] = ynode[iper-1+i];
    xsup[3] = xnode[nreg] + (xnode[1] - xnode[0]);
    ysup[3] = ynode[1];
  } else {
    for (int i = 0; i < 4; ++i) xsup[i] = xnode[iper-1+i];
    for (int i = 0; i < 4; ++i) ysup[i] = ynode[iper-1+i];
  }

  return eval_lagrange_poly(4, xsup, ysup, xper);
}

// Very straightforward, inefficient impl.
void periodic_cubic1d_make_translation_matrix (
  const Int nreg, const Real* const xnode, const Real& xoffset, Real* op, Real* wrk)
{
  for (Int i = 0; i < nreg; ++i) wrk[i] = 0;
  for (Int si = 0; si < nreg; ++si) {
    wrk[si] = 1;
    if (si == 0) wrk[nreg] = 1;
    for (Int ti = 0; ti < nreg; ++ti)
      op[nreg*si + ti] = periodic_cubic1d_interp(nreg, xnode, wrk, xnode[ti] + xoffset);
    wrk[si] = 0;
    if (si == 0) wrk[nreg] = 0;
  }
}

// Demo a configuration that has associated max eig amp >= 1 + 1e-3.
static Real cubic1d_demo_unstable_problem () {
  const Real x[] = {0, 0.11242, 0.44817, 0.78392, 0.88737, 1};
  const Int nreg = sizeof(x)/sizeof(*x) - 1;
  const Real xoffset = 0.33575;
  std::vector<Real> op(nreg*nreg), wrk(nreg+1);
  periodic_cubic1d_make_translation_matrix(nreg, x, xoffset, op.data(), wrk.data());
  const auto mea = calc_max_eig_amp(nreg, op.data(), nreg, wrk);
  return mea;
}

static void cubic1d_unittest () {
  const auto eps = std::numeric_limits<Real>::epsilon();

  {
    require(equal(calc_periodic_x(-0.1, 1.3, 0.7), 0.7));
    require(equal(calc_periodic_x(-0.1, 1.3, 1.3), 1.3));
    require(equal(calc_periodic_x(-0.1, 1.3, -0.1), -0.1));
    require(almost_equal(calc_periodic_x(-0.1, 1.3, 1.4), 0, eps));
    const auto x = calc_periodic_x(-0.1, 1.3, 2.7);
    require(almost_equal(x, -0.1, eps) || almost_equal(x, 1.2, eps));
    require(almost_equal(calc_periodic_x(1.1, 1.3, 1.4), 1.2, eps));
  }

  { // In the interior, a cubic is recovered exactly.
    const auto f = [&] (const Real& x) { return ((((-1.2*x) + 0.1)*x) - 0.3)*x + 11; };
    const Int nreg = 5, nsamp = 100;
    std::vector<Real> x(nreg+1), y(nreg+1);
    for (Int i = 0; i <= nreg; ++i) x[i] = 2*urand() - 1;
    std::sort(x.begin(), x.end());
    for (Int i = 0; i < nreg; ++i) y[i] = f(x[i]);
    y[nreg] = y[0]; // y[nreg] doesn't influence this test.
    const auto D = x[nreg-2] - x[1];
    for (Int i = 0; i < nsamp; ++i) {
      const auto xs = x[1] + i*D/(nsamp - 1);
      const auto ys = periodic_cubic1d_interp(nreg, x.data(), y.data(), xs);
      const auto yt = f(xs);
      require(almost_equal(ys, yt, 10*eps));
    }
  }

  { // At each sample point, error decreases at least as fast as nreg^-4 in the
    // 1-norm, nreg^-3 pointwise.
    const Real xl = -1.2, xu = 4.2, L = xu - xl;
    const auto f = [&] (const Real& x) { return std::cos(2*M_PI*(x - xl)/L); };
    const Int nsamp = 1111;
    std::vector<Real> xs(nsamp);
    std::vector<std::vector<Real> > yss(2);
    for (Int i = 0; i < nsamp; ++i) xs[i] = xl + i*L/(nsamp - 1);
    for (Int nreg : {31, 80, 211}) {
      Int cnt = 0;
      for (int refine = 0; refine < 2; ++refine) {
        nreg *= 2;
        std::vector<Real> x(nreg+1), y(nreg+1);
        for (Int i = 0; i <= nreg; ++i) x[i] = xl + i*L/nreg;
        for (Int i = 0; i <= nreg; ++i) y[i] = f(x[i]);
        auto& ys = yss[refine];
        ys.resize(nsamp);
        for (Int i = 0; i < nsamp; ++i)
          ys[i] = periodic_cubic1d_interp(nreg, x.data(), y.data(), xs[i]);
      }
      Real err[2] = {0};
      for (Int i = 0; i < nsamp; ++i) {
        const auto yt = f(xs[i]);
        const auto e1 = std::abs(yss[0][i] - yt);
        const auto e2 = std::abs(yss[1][i] - yt);
        require(e1 >= e2);
        if ( ! (e1 == e2 || e1 >= 8*e2))
          ++cnt;
        if (i < nsamp-1) {
          err[0] += e1;
          err[1] += e2;
        }
      }
      // A few points not showing 3rd-order convergence is OK, since this means
      // the less-accurate solution is more accurate than anticipated at these
      // points.
      require(cnt < 0.1*nsamp);
      require(err[0] > 15.9*err[1]);
      require(err[0] < 16.1*err[1]);
    }
  }

  { // The matrix gives the same answer as calls to periodic_cubic1d_interp.
    const Int nreg = 17;
    std::vector<Real> x(nreg+1), op(nreg*nreg), wrk(nreg+1),
      ys(nreg+1), yt1(nreg+1), yt2(nreg+1);

    for (Int i = 0; i <= nreg; ++i) x[i] = 2*urand() - 1;
    std::sort(x.begin(), x.end());
    for (Int i = 0; i <  nreg; ++i) ys[i] = 2*urand() - 1;
    ys[nreg] = ys[0];
    
    for (const auto xoffset : {0.0, 0.01, -0.02, 0.1, -0.42, 1.7, -4.2}) {
      periodic_cubic1d_make_translation_matrix(nreg, x.data(), xoffset, op.data(), wrk.data());
      for (Int i = 0; i < nreg; ++i) yt1[i] = 0;
      for (Int j = 0; j < nreg; ++j)
        for (Int i = 0; i < nreg; ++i)
          yt1[i] += op[nreg*j + i]*ys[j];
      for (Int i = 0; i < nreg; ++i)
        yt2[i] = periodic_cubic1d_interp(nreg, x.data(), ys.data(), x[i] + xoffset);
      for (Int i = 0; i < nreg; ++i)
        require(almost_equal(yt1[i], yt2[i], 10*eps));
    }
  }

  { // The matrix for a uniform grid has max eig 1.
    const Int nreg = 6;
    std::vector<Real> x(nreg+1), op(nreg*nreg), wrk(nreg+1);
    for (Int i = 0; i <= nreg; ++i) x[i] = i;
    const Real xoffset = 1.2;
    periodic_cubic1d_make_translation_matrix(nreg, x.data(), xoffset, op.data(), wrk.data());
    const auto mea = calc_max_eig_amp(nreg, op.data(), nreg, wrk);
    require(almost_equal(mea, 1, 2*eps));
  }

  { // Cubic ISL on the uniform-grid periodic translation problem has OOA 3.
    const Real xl = -4.2, xu = 1.7, L = xu - xl;
    const auto f = [&] (const Real& x) { return std::cos(2*M_PI*(x - xl)/L); };
    const auto error = [&] (const std::vector<Real>& x, const std::vector<Real>& y) {
      Real num = 0, den = 0;
      // Ignore the last point, which is periodic; we don't update it.
      for (size_t i = 0; i < y.size(); ++i) {
        const auto yt = f(x[i]);
        num += square(y[i] - yt);
        den += square(yt);
      }
      return std::sqrt(num/den);
    };
    Int nstep = 37;
    Int nreg = 20;
    Real e[2];
    for (Int refine = 0; refine < 2; ++refine) {
      nreg *= 2;
      nstep *= 2;
      const Real xoffset = L/nstep;
      std::vector<Real> x(nreg+1), op(nreg*nreg), wrk(nreg+1);
      std::vector<Real> ys[2];
      // Size is nreg rather than nreg+1 because we don't maintain the periodic
      // point y[nreg] = y[0].
      for (Int k = 0; k < 2; ++k) ys[k].resize(nreg);
      // Make the space-time operator A.
      for (Int i = 0; i <= nreg; ++i) x[i] = xl + i*L/nreg;
      periodic_cubic1d_make_translation_matrix(nreg, x.data(), xoffset, op.data(), wrk.data());
      // Initial condition.
      Int i0 = 0, i1 = 1;
      for (Int i = 0; i < nreg; ++i) ys[i0][i] = f(x[i]);
      for (Int si = 0; si < nstep; ++si) {
        // y1 = A y0
        for (Int i = 0; i < nreg; ++i) ys[i1][i] = 0;
        for (Int j = 0; j < nreg; ++j)
          for (Int i = 0; i < nreg; ++i)
            ys[i1][i] += op[nreg*j + i]*ys[i0][j];
        // At a time t < T, check that the error is large.
        if (si == nstep/2) require(error(x, ys[i1]) > 1);
        std::swap(i0, i1);
      }
      e[refine] = error(x, ys[i0]);
    }
    // At time T, the error decreases with OOA 3.
    require(e[0] > 7.95*e[1]);
    require(e[0] < 8.05*e[1]);
    // Check that we solved the problem to reasonable accuracy.
    require(e[1] <= 1e-3);
  }

  // The primary purpose of periodic_cubic1d_make_translation_matrix is to
  // demonstrate that there is a 1D periodic translation problem for which the
  // associated classical cubic ISL space-time matrix has maximum eigenvalue
  // amplitude > 1. In this demo, it is > 1 + 1e-3.
  require(cubic1d_demo_unstable_problem() >= 1 + 1e-3);
}

struct NonUniMesh1d {
  using Array = std::vector<Real>;

  NonUniMesh1d (const Real* xb, const Int ne)
    : ne_(ne), L_(xb[ne] - xb[0]), xb_(ne+1)
  {
    std::copy(xb, xb+ne+1, xb_.begin());
  }

  Int get_ne () const { return ne_; }

  const Array& get_xb() const { return xb_; }

  Real to_periodic (const Real& x) const {
    if (x >= xb_[0] && x <= xb_[ne_]) return x;
    auto y = (x - xb_[0])/L_;
    y = y - std::floor(y);
    return xb_[0] + y*L_;
  }

  Int in_cell (const Real& x) const {
    const auto xp = to_periodic(x);
    if (xp <= xb_[0]) return 0;
    const Int iper = std::lower_bound(xb_.begin(), xb_.end(), xp) - xb_.begin() - 1;
    assert(iper >= 0  );
    assert(iper <  ne_);
    assert(xp >= xb_[iper  ]);
    assert(xp <= xb_[iper+1]);
    return iper;
  }

  Real to_physical (const Int& ie) const { return xb_[ie]; }

  static Int unittest () {
    const Real eps = std::numeric_limits<Real>::epsilon();
    Int ne = 0;
    std::vector<Real> xb({-1, 1, 5.});
    NonUniMesh1d m(xb.data(), 2);
    if (m.in_cell(1.1) != 1) ++ne;
    if (m.in_cell(-1) != 0) ++ne;
    if (reldif(m.to_periodic(5.1), -0.9) > 10*eps) ++ne;
    if (m.in_cell(5.1) != 0) ++ne;
    if (m.in_cell(-1.1) != 1) ++ne;
    return ne;
  }

private:
  const Int ne_;
  const Real L_;
  Array xb_;
};

class NonUniMesh2d {
  const NonUniMesh1d mx_, my_;
  const Int nx_;

public:
  NonUniMesh2d (const Real* xb, const Int nx, const Real* yb, const Int ny)
    : mx_(xb, nx), my_(yb, ny), nx_(nx)
  {}

  Int get_ne () const { return mx_.get_ne() * my_.get_ne(); }

  const NonUniMesh1d& get_mx () const { return mx_; }
  const NonUniMesh1d& get_my () const { return my_; }

  void to_periodic (const Real& x, const Real& y,
                    Real& xp, Real& yp) const {
    xp = mx_.to_periodic(x);
    yp = my_.to_periodic(y);
  }

  void ncell (Int& nx, Int& ny) const {
    nx = nx_;
    ny = my_.get_ne();
  }

  Int in_cell (const Real& x, const Real& y) const {
    const Int ix = mx_.in_cell(x), iy = my_.in_cell(y);
    return iy*nx_ + ix;
  }

  void to_physical (const Int& ci, Real& x, Real& y) const {
    const Int yci = ci / nx_;
    const Int xci = ci % nx_;
    x = mx_.to_physical(xci);
    y = my_.to_physical(yci);
  }

  static Int unittest () {
    const Real eps = std::numeric_limits<Real>::epsilon();
    Int ne = 0;
    std::vector<Real> xb({-1, 1, 5.}), yb({-1, 1, 6.});
    NonUniMesh2d m(xb.data(), 2, yb.data(), 2);
    if (m.in_cell(1.1, 1.1) != 3) ++ne;
    if (m.in_cell(5, 6) != 3) ++ne;
    if (m.in_cell(-1, 1) != 0) ++ne;
    if (m.in_cell(1, 1) != 0) ++ne;
    Real xp, yp;
    m.to_periodic(5.1, -1.2, xp, yp);
    if (reldif(xp, -0.9) > 10*eps) ++ne;
    if (reldif(yp, 5.8) > 10*eps) ++ne;
    if (m.in_cell(5.1, 6.1) != 0) ++ne;
    if (m.in_cell(-1.1, -1.1) != 3) ++ne;
    return ne;
  }
};

// Sparse matrix data structures and operations.
struct SparseTriple {
  Int m, n;
  std::vector<Int> rp, ci; // row pointer, column index
  std::vector<Real> d; // matrix entries
};

// d is row major.
static void sparse2dense (const SparseTriple& s, std::vector<Real>& d) {
  d.resize(s.m*s.n, 0);
  for (Int r = 0; r < s.m; ++r)
    for (Int j = s.rp[r]; j < s.rp[r+1]; ++j)
      d[s.n*r + s.ci[j]] = s.d[j];
}

static void apply (const SparseTriple& s, const std::vector<Real>& x,
                   std::vector<Real>& y) {
  assert(s.n == static_cast<Int>(x.size()));
  y.resize(s.m);
  for (Int r = 0; r < s.m; ++r) {
    Real yr = 0;
    for (Int j = s.rp[r]; j < s.rp[r+1]; ++j)
      yr += s.d[j]*x[s.ci[j]];
    y[r] = yr;
  }
}

// Fill the polynomial interpolant's support with periodically unwrapped
// coordinate values.
static void fill_lag_coord (const NonUniMesh1d::Array& xb, const Int& np,
                            const Int& cell, Real* coord) {
  const Int nx = xb.size();
  const Int ne = nx - 1;
  const Real L = xb[ne] - xb[0];
  for (Int i = 0; i < np; ++i) {
    const Int k = cell + i;
    if (k < 0)
      coord[i] = xb[(k + ne) % ne] - L;
    else if (k >= ne)
      coord[i] = xb[k % ne] + L;
    else
      coord[i] = xb[k];
  }
#ifndef NDEBUG
  for (Int i = 1; i < np; ++i) assert(coord[i] > coord[i-1]);
#endif
}

// Order the degrees of freedom in the operator.
static Int dof (const Int& nx, const Int& ny,
                const Int& xe, const Int& ye) {
  return (ye % ny)*nx + (xe % nx);
}

template <typename Function>
void make_ccsl_op_nondiv2d (const NonUniMesh2d& mesh, const Int np,
                            const Function& integrate,
                            const Real& dt, SparseTriple& s)
{
  const Int os = (np-1)/2;
  Int nex, ney;
  mesh.ncell(nex, ney);
  const Int ne = nex*ney, n = ne;
  const auto& xb = mesh.get_mx().get_xb();
  const auto& yb = mesh.get_my().get_xb();

  s.m = s.n = ne;
  s.rp.resize(ne+1, 0);

  for (Int tie = 0; tie < ne; ++tie) {
    const Int tye = tie / nex, txe = tie % nex;
    const Int tdof = dof(nex, ney, txe, tye);
    assert(tdof >= 0 && tdof < n);

    Real x0, y0;
    mesh.to_physical(tie, x0, y0);
    Real tx, ty;
    integrate(x0, y0, dt, tx, ty);

    const Int sie = mesh.in_cell(tx, ty);
    const Int sxe0 = sie % nex, sye0 = sie / nex;

    Real txp, typ;
    mesh.to_periodic(tx, ty, txp, typ);
    Real xv[12], yv[12], lag_coord[12];
    assert(np <= 12);
    fill_lag_coord(xb, np, sxe0 - os, lag_coord);
    assert(txp >= lag_coord[0] && txp <= lag_coord[np-1]);
    eval_lagrange_poly_basis(np, lag_coord, txp, xv);
    fill_lag_coord(yb, np, sye0 - os, lag_coord);
    assert(typ >= lag_coord[0] && typ <= lag_coord[np-1]);
    eval_lagrange_poly_basis(np, lag_coord, typ, yv);

    Int k = s.rp[tie];
    for (Int dy = -os; dy <= os+1; ++dy) {
      const Int sye = (sye0 + dy + ney) % ney;
      for (Int dx = -os; dx <= os+1; ++dx) {
        const Int sxe = (sxe0 + dx + nex) % nex;
        const Int sdof = dof(nex, ney, sxe, sye);
        assert(sdof >= 0 && sdof < n);
        s.ci.push_back(sdof);
        s.d.push_back(xv[dx+os]*yv[dy+os]);
        k++;
      }
    }
    s.rp[tie+1] = k;
  }

  assert(s.d.size() == s.ci.size());
  assert(s.rp[s.m] == static_cast<Int>(s.d.size()));
}

// Check that the 'integrate' function provides doubly-periodic outputs over the
// domain of the mesh.
template <typename Function>
void check_periodicity (const NonUniMesh2d& mesh, const Real& dt, const Function& integrate) {
  static const Real rtol = 1e2*std::numeric_limits<Real>::epsilon();
  static const Real atol = std::numeric_limits<Real>::epsilon();

  Int nex, ney;
  mesh.ncell(nex, ney);
  const auto& xb = mesh.get_mx().get_xb();
  const auto& yb = mesh.get_my().get_xb();

  // Check that the domain is [0,1]^2.
  require(std::abs(xb[nex] - xb[0] - 1) <= atol);
  require(std::abs(yb[ney] - yb[0] - 1) <= atol);

  // Relative/absolute error checks.
  const auto check_error = [=] (Real x0, Real tx0, Real x1, Real tx1) {
    require(reldif(tx0 - x0, tx1 - x1) <= rtol ||
            std::abs((tx0 - x0) - (tx1 - x1)) <= atol);
  };

  // Run over the sides of the domain and check periodicity.
  const auto check = [&] (Real x0, Real x1, Real y0, Real y1) {
    Real tx0, ty0, tx1, ty1;
    integrate(x0, y0, dt, tx0, ty0);
    integrate(x1, y1, dt, tx1, ty1);
    check_error(x0, tx0, x1, tx1);
    check_error(y0, ty0, y1, ty1);
  };
  Int n = 7*std::max(nex, ney);
  for (int i = 0; i <= n; ++i) {
    const Real a = Real(i)/n, x = a*xb[0] + (1-a)*xb[nex];
    check(x, x, yb[0], yb[ney]);
  }
  for (int i = 0; i <= n; ++i) {
    const Real a = Real(i)/n, y = a*yb[0] + (1-a)*yb[ney];
    check(xb[0], xb[nex], y, y);
  }
}

static void setup_demo_problem (const std::vector<Real>& xb, const std::vector<Real>& yb,
                                const Real dt, SparseTriple& s, const Int np,
                                const bool check = false) {
  const Int nex = xb.size() - 1, ney = yb.size() - 1;
  // Shear (nondivergent) flow. Parameter values were obtained from a search for
  // an unstable operator.
  const auto integrate = [&] (const Real& x0, const Real& y0, const Real& dt,
                              Real& xf, Real& yf) {
    const auto speed = 1 + std::cos(2*M_PI*(0.342 + x0 - y0));
    xf = x0 + speed*dt;
    yf = y0 + speed*dt;
  };
  NonUniMesh2d mesh(xb.data(), xb.size()-1, yb.data(), yb.size()-1);
  if (check) check_periodicity(mesh, dt, integrate);
  make_ccsl_op_nondiv2d(mesh, np, integrate, dt, s);
}

static Real cubic2d_demo_unstable_problem (const bool unstable = true) {
  const Int nex = 15, ney = 13, ne = nex*ney;
  const Real dt = unstable ? 0.2761 : 0.1;
  std::vector<Real> xb(nex+1), yb(ney+1), op(ne*ne), wrk;
  const Real dx = 1.0/nex, dy = 1.0/ney;
  for (Int i = 0; i <= nex; ++i) xb[i] = i*dx;
  for (Int i = 0; i <= ney; ++i) yb[i] = i*dy;
  SparseTriple s;
  setup_demo_problem(xb, yb, dt, s, 4, true);
  sparse2dense(s, op);
  const auto mea = calc_max_eig_amp(ne, op.data(), ne, wrk);
  return mea;
}

static void set_ic (const std::vector<Real>& xb, const std::vector<Real>& yb,
                    std::vector<Real>& z) {
  const Int nex = xb.size() - 1, ney = yb.size() - 1, ne = nex*ney;
  z.resize(ne);
  for (Int iy = 0; iy < ney; ++iy) {
    const Real fy = std::cos(2*M_PI*(yb[iy] - yb[0])/(yb[ney] - yb[0]));
    for (Int ix = 0; ix < nex; ++ix)
      z[nex*iy + ix] = fy*std::cos(2*M_PI*(xb[ix] - xb[0])/(xb[nex] - xb[0]));
  }
}

static Real calc_l2_reldif (const std::vector<Real>& a, const std::vector<Real>& b) {
  Real num = 0, den = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    num += square(a[i] - b[i]);
    den += square(a[i]);
  }
  return std::sqrt(num/den);
}

// The example problem generically has order of accuracy np-1 as long as shear
// is not permitted to continue for too long. Demonstrate this:
//   If reverse, integrate a smooth field forward for time T/2, then backward to
// time T, and compare the result with the initial condition.
//   If not reverse, integrate forward for time T, then compare with the
// analytical solution at time T.
static Real measure_ooa (const Int np, const bool reverse) {
  Int nex = 35, ney = 32, nstep = 22;
  const Real x0 = -0.4, y0 = 0.2;
  Real dt = 0.012, err[2];
  const Real T = dt*nstep;
  for (Int refine = 0; refine < 2; ++refine) {
    // Refinement parameters.
    nex *= 2;
    ney *= 2;
    nstep *= 2;
    dt /= 2;
    const Int ne = nex*ney;
    // Make forward and backward ops.
    std::vector<Real> xb(nex+1), yb(ney+1);
    for (Int i = 0; i <= nex; ++i) xb[i] = x0 + Real(i)/nex;
    for (Int i = 0; i <= ney; ++i) yb[i] = y0 + Real(i)/ney;
    SparseTriple op1, op2;
    setup_demo_problem(xb, yb, dt, op1, np, refine == 0);
    setup_demo_problem(xb, yb, reverse ? -dt : T, op2, np, refine == 0);
    // Initial conditions.
    std::vector<Real> z0, zs[2];
    set_ic(xb, yb, z0);
    zs[0] = z0;
    zs[1].resize(z0.size());
    // Time step.
    Int i0 = 0, i1 = 1;
    for (Int ti = 0; ti < nstep; ++ti) {
      apply(op1, zs[i0], zs[i1]);
      std::swap(i0, i1);
    }
    if (reverse) {
      for (Int ti = 0; ti < nstep; ++ti) {
        apply(op2, zs[i0], zs[i1]);
        std::swap(i0, i1);
      }
      err[refine] = calc_l2_reldif(z0, zs[i0]);
    } else {
      apply(op2, z0, zs[i1]);
      err[refine] = calc_l2_reldif(zs[i1], zs[i0]);
    }
  }
  return std::log2(err[0]/err[1]);
}

void cubic2d_unittest () {
  const auto eps = std::numeric_limits<Real>::epsilon();
  require(NonUniMesh1d::unittest() == 0);
  require(NonUniMesh2d::unittest() == 0);
  // The primary purpose of make_ccsl_op_nondiv2d is to demonstrate that there
  // is a 2D periodic nondivergent-flow problem, with uniform grid, for which
  // the associated classical cubic ISL space-time matrix has maximum eigenvalue
  // amplitude > 1. In this demo, it is > 1 + 1e-2.
  require(cubic2d_demo_unstable_problem() >= 1 + 1e-2);
  // Make sure it's stable if, for example, dt is not from the seach for
  // unstable parameter values.
  require(almost_equal(cubic2d_demo_unstable_problem(false), 1, 50*eps));
  // Order of accuracy is np-1.
  for (const Int np : {4, 6, 8})
    for (const bool reverse : {true, false}) {
      const auto ooa = measure_ooa(np, reverse);
      const Real d = 0.025; // permit 2.5% deviation from theoretical OOA
      require(ooa > (1 - d)*(np - 1));
      require(ooa < (1 + d)*(np - 1));
    }
}

int main (int argc, char** argv) {
  // Show a 1D periodic translation problem on a nonuniform grid for which the
  // associated classical cubic ISL space-time matrix has maximum eigenvalue
  // amplitude > 1 + 1e-3.
  cubic1d_unittest();
  // Show a 2D periodic nondivergent-flow problem on a uniform grid for which
  // the associated classical cubic ISL space-time matrix has maximum eigenvalue
  // amplitude > 1 + 1e-2.
  cubic2d_unittest();
}
