#include "slmmir_physgrid.hpp"

#include "slmm_defs.hpp"
#include "slmm_mesh.hpp"
#include "slmm_spf.hpp"
#include "slmm_gll.hpp"
#include "slmm_io.hpp"
#include "slmm_nla.hpp"
#include "slmm_debug.hpp"
#include "slmm_basis_reduced.hpp"

#include "slmmir_util.hpp"

namespace pg {

// Cell-local limiter. N.B. We don't handle rho <= 0 in this study.
void limit (
  const Limiter::Enum limiter, const Real qlo, const Real qhi,
  const Array<Real>& twt, const Real* tmetdet, Real* trho, Real* tq)
{
  if (limiter == Limiter::none) return;
  // The particular limiter is unimportant in this study, so just use CAAS.
  if (limiter != Limiter::caas)
    throw std::runtime_error("Only CAAS is supported");
  const Int nt = twt.size();
  Real a[GLL::np_max*GLL::np_max];
  Real Qmass = 0;
  for (Int i = 0; i < nt; ++i) {
    a[i] = twt[i]*tmetdet[i]*trho[i];
    Qmass += a[i]*tq[i];
  }
  Real qstar[GLL::np_max*GLL::np_max];
  for (Int i = 0; i < nt; ++i)
    qstar[i] = tq[i];
  spf::clip_and_sum(twt.size(), nullptr, a, Qmass, &qlo, &qhi, true, qstar, tq);
}

Fv2Gll::Type::Enum Fv2Gll::Type::convert (const std::string& s) {
  if (s == "idem") return Type::idem;
  if (s == "l2") return Type::l2;
  if (s == "l2ep") return Type::l2ep;
  if (s == "elrecon") return Type::elrecon;
  SIQK_THROW_IF(true, "Fv2Gll::convert: not impled for " << s);
}

std::string Fv2Gll::Type::convert (const Fv2Gll::Type::Enum e) {
  switch (e) {
  case Type::idem: return "idem";
  case Type::l2: return "l2";
  case Type::l2ep: return "l2ep";
  case Type::elrecon: return "elrecon";
  default: SIQK_THROW_IF(true, "Fv2Gll::convert: not impled for " << e);
  }
}

Fv2Gll::Ptr Fv2Gll::create (Fv2Gll::Type::Enum type, Int np, Int nphys,
                            const Basis::ConstPtr& b) {
  if (nphys < 3 && type == Fv2Gll::Type::elrecon) {
    printf("WARNING: Switching from elrecon to idem since nf < 3.\n");
    type = Fv2Gll::Type::idem;
  }
  switch (type) {
  case Fv2Gll::Type::idem: return std::make_shared<IdemFv2Gll>(np, nphys, b);
  case Fv2Gll::Type::l2: return std::make_shared<L2Fv2Gll>(np, nphys, b);
  case Fv2Gll::Type::l2ep: return std::make_shared<L2ExceptPerimFv2Gll>(np, nphys, b);
  case Fv2Gll::Type::elrecon: return std::make_shared<ElemLclReconFv2Gll>(np, nphys, b);
  default: SIQK_THROW_IF(true, "Fv2Gll::create: not impled for " << type);
  }
}

// 3 works but requires a limiter to be applied to each region (separately) in
// reconstruct_nphys1. 2 gives approx the same OOA as 3 but requires no
// limiter. So let's go with that for now.
const Int IdemFv2Gll::np_nphys1 = 2;

void Remap::init (Int np_, Int nphys_, const Basis::ConstPtr& b_) {
  np = np_;
  nphys = nphys_;
  b = b_;

  w_dd.optclear_and_resize(np*np);

  const Real* wt;
  b->get_w(np, wt);

  // Diagonal (lumped) mass matrix. In all the integrals, I ignore the 0.25
  // factor as it cancels out in the projections.
  for (Int i = 0; i < np; ++i)
    for (Int j = 0; j < np; ++j)
      w_dd[np*i + j] = wt[i]*wt[j];
}

// Full mixed mass matrix, column-major. Can be impl'ed with Homme's
// subcell_integration.
void Remap::init_M_dp (const Int np, const Int nphys, const Basis& b,
                       Array<Real>& M_dp) {
  const Int nphys2 = nphys*nphys, np2 = np*np;
  M_dp.optclear_and_resize(np2*nphys2, 0);
  Basis::compute_integrals_over_subcells_2d(b, np, nphys, M_dp.data());
}

void Gll2Fv::init (Int np, Int nphys, const Basis::ConstPtr& basis) {
  Remap::init(np, nphys, basis);
  init_matrices();
}

void Gll2Fv::init_matrices () {
  const Int nphys2 = nphys*nphys;

  M_pp.optclear_and_resize(nphys2, 0);

  // Diagonal FV mass matrix. Note the sum is 4 since the ref elem is [-1,1]^2
  // and we don't multiply integrals by 0.25.
  for (Int i = 0; i < nphys2; ++i)
    M_pp[i] = square(2.0/nphys);

  init_M_dp(np, nphys, *b, M_dp);
}

void Gll2Fv::remapd (const Real* gll_metdet, const Real* fv_metdet,
                     const Real* d, Real* p) const {
  const Int nphys2 = nphys*nphys, np2 = np*np;
  for (Int pi = 0; pi < nphys2; ++pi) {
    Real accum = 0;
    for (Int di = 0; di < np2; ++di)
      accum += M_dp[np2*pi + di]*d[di]*gll_metdet[di];
    p[pi] = accum/(M_pp[pi]*fv_metdet[pi]);
  }  
}

void Gll2Fv::remap (const Real* gll_metdet, const Real* fv_metdet,
                    const Limiter::Enum limiter,
                    const Real* drho, const Real* dq,
                    Real* prho, Real* pq, const bool remap_rho) const {
  const Int nphys2 = nphys*nphys, np2 = np*np;
  if (remap_rho) remapd(gll_metdet, fv_metdet, drho, prho);
  Real wrk[GLL::np_max*GLL::np_max];
  for (Int i = 0; i < np2; ++i)
    wrk[i] = drho[i]*dq[i];
  remapd(gll_metdet, fv_metdet, wrk, pq);
  for (Int i = 0; i < nphys2; ++i)
    pq[i] /= prho[i];
  // In the GLL -> FV direction, limiter bounds come from just this element.
  Real qlo = dq[0], qhi = dq[0];
  for (Int i = 1; i < np2; ++i) qlo = std::min(qlo, dq[i]);
  for (Int i = 1; i < np2; ++i) qhi = std::max(qhi, dq[i]);
  limit(limiter, qlo, qhi,
        M_pp, fv_metdet, prho, pq);
}

void Fv2Gll
::remap (const Real* gll_metdet, const Real* fv_metdet,
         // In the FV -> GLL direction, we need to expand the bounds
         // to those coming from neighboring cells.
         const Limiter::Enum limiter, const Real pqlo, const Real pqhi,
         const Real* prho, const Real* pq,
         Real* drho, Real* dq, const bool remap_rho) const {
  const Int nphys2 = nphys*nphys, np2 = np*np;
  if (remap_rho) remapd(gll_metdet, fv_metdet, prho, drho);
  Real wrk[GLL::np_max*GLL::np_max];
  for (Int i = 0; i < nphys2; ++i)
    wrk[i] = prho[i]*pq[i];
  remapd(gll_metdet, fv_metdet, wrk, dq);
  for (Int i = 0; i < np2; ++i)
    dq[i] /= drho[i];
  if (nphys > 1)
    limit(limiter, pqlo, pqhi,
          w_dd, gll_metdet, drho, dq);
}

void IdemFv2Gll::init (Int np, Int nphys, const Basis::ConstPtr& basis) {
  Remap::init(np, nphys, basis);
  // We decided to use just nphys.
  npi = nphys == 1 ? np_nphys1 : nphys; // std::min(np, std::max(3, nphys));
  init_matrices();
}

static void apply_idem_op (const Int np, const Int npi, const Int nphys,
                           const Array<Real>& R, const Array<Real>& tau,
                           const Array<Real>& M_dp, const Array<Real>& npi_to_np,
                           const Array<Real>& M_dd, const Array<Real>& M_pp,
                           const Real* const p, Real* const d) {  
  const Int nphys2 = nphys*nphys, npi2 = npi*npi, np2 = np*np;
  assert(npi <= np);
  assert(nphys2 <= npi2);
  // Constrained projection, as described in init_matrices.
  Real wrk[GLL::np_max*GLL::np_max];
  //   wrk = M_pp (fv_metdet p)
  for (Int pi = 0; pi < nphys2; ++pi)
    wrk[pi] = M_pp[pi]*p[pi];
  if (nphys < npi) {
    //   Solve S wrk = wrk using S = R'R.
    dtrtrs('u', 't', 'n', nphys2, 1, R.data(), npi2, wrk, nphys2);
    dtrtrs('u', 'n', 'n', nphys2, 1, R.data(), npi2, wrk, nphys2);
    //   d = M_dp wrk
    for (Int di = 0; di < npi2; ++di)
      d[di] = 0;
    for (Int pi = 0; pi < nphys2; ++pi)
      for (Int di = 0; di < npi2; ++di)
        d[di] += M_dp[npi2*pi + di]*wrk[pi];
  } else {
    Real wrk2[GLL::np_max*GLL::np_max];
    dtrsm('l', 'u', 't', 'n', nphys2, 1, 1.0, R.data(), npi2, wrk, nphys2);
    dormqr('l', 'n', nphys2, 1, nphys2, R.data(), npi2, tau.data(), wrk, nphys2,
           wrk2, np2);
    for (Int di = 0; di < npi2; ++di)
      d[di] = wrk[di];
  }
  if (npi < np) {
    if (nphys == npi) {
      for (Int di = 0; di < npi2; ++di)
        wrk[di] = d[di];
    } else {
      //   wrk = inv(M_dd) d
      for (Int di = 0; di < npi2; ++di)
        wrk[di] = d[di]/M_dd[di];
    }
    // Interpolate from npi to np; if npi == np, this is just the Id matrix.
    for (Int i = 0; i < np2; ++i) {
      Real accum = 0;
      for (Int ii = 0; ii < npi2; ++ii)
        accum += npi_to_np[npi2*i + ii]*wrk[ii];
      d[i] = accum;
    }
  } else if (nphys < npi) {
    // Do d = inv(M_dd) d and divide out the ref -> sphere jacobian.
    for (Int di = 0; di < np2; ++di)
      d[di] /= M_dd[di];
  }
}

void IdemFv2Gll::init_matrices () {
  const Int nphys2 = nphys*nphys, npi2 = npi*npi, np2 = np*np;
  Array<Real> M_pp, M_dd, M_dp, R, tau;

  M_pp.optclear_and_resize(nphys2, 0);

  // diagonal FV mass matrix
  for (Int i = 0; i < nphys2; ++i)
    M_pp[i] = square(2.0/nphys);

  GLL g_npi; // GLL here and not basis b/c this is the intermediate basis

  if (nphys < npi) {
    assert(0);
    M_dd.optclear_and_resize(npi2, 0);
    const Real* wt;
    g_npi.get_w(npi, wt);
    // diagonal (lumped) GLL mass matrix
    for (Int i = 0; i < npi; ++i)
      for (Int j = 0; j < npi; ++j)
        M_dd[npi*i + j] = wt[i]*wt[j];
  }

  init_M_dp(npi, nphys, g_npi, M_dp);

  if (false && npi == nphys) {
    const Int nwrk = 10*npi2*npi2;
    Array<Real> A(npi2*npi2), s(npi2), wrk(nwrk);
    for (Int i = 0; i < npi2*npi2; ++i) A[i] = M_dp[i];
    dgesvd('n', 'n', npi2, npi2, A.data(), npi2, s.data(),
           nullptr, 1, nullptr, 1,
           wrk.data(), nwrk);
    printf("cond %1.3e\n", s[0]/s[npi2-1]);
  }

  /* We want to solve
         min_d 1/2 d'M_dd d - d' M_dp p
          st   M_dp' d = M_pp p,
     which gives
         [M_dd -M_dp] [d] = [M_dp p]
         [M_dp'  0  ] [y]   [M_pp p].
     Recall M_dd, M_pp are diag. Let
         S = M_dp' inv(M_dd) M_dp.
     Then
         d = inv(M_dd) M_dp inv(S) M_pp p.
     Compute the QR factorization sqrt(inv(M_dd)) M_dp = Q R so that S = R'R.
       If nphys = npi, then the problem reduces to
         M_dp' d = M_pp p.
     M_dp is symmetric. Compute
         M_dp = Q R
     and later solve
         R'Q' d = M_pp p.
   */
  R.optclear_and_resize(npi2*nphys2);
  for (Int pi = 0; pi < nphys2; ++pi)
    for (Int di = 0; di < npi2; ++di)
      R[npi2*pi + di] = (nphys < npi ?
                         M_dp[npi2*pi + di]/std::sqrt(M_dd[di]) :
                         M_dp[npi2*pi + di]);
  tau.optclear_and_resize(npi2*nphys2);
  Array<Real> wrk(npi2*nphys2);
  dgeqrf(npi2, nphys2, R.data(), npi2, tau.data(), wrk.data(), npi2*nphys2);

  // interpolation matrix for npi -> np
  if (b->max_degree(np) < npi-1)
    printf("WARNING: npi_to_np matrix is not mass conserving: %d vs %d.\n",
           b->max_degree(np), npi-1);
  build_npi_to_np_matrix(npi, np, *b, npi_to_np);

  op_p_to_d.optclear_and_resize(np2*nphys2);
  Array<Real> e(nphys2, 0), col(np2);
  for (Int i = 0; i < nphys2; ++i) {
    e[i] = 1;
    apply_idem_op(np, npi, nphys, R, tau, M_dp, npi_to_np, M_dd, M_pp,
                  e.data(), col.data());
    for (Int j = 0; j < np2; ++j) op_p_to_d[nphys2*j + i] = col[j];
    e[i] = 0;
  }
}

static void build_interp_matrix (const Int np_from, const Int np_to,
                                 const Basis& b, Array<Real>& interp) {
  const Int npf2 = np_from*np_from, npt2 = np_to*np_to;
  const Real* fx, * tx;
  b.get_x(np_from, fx);
  b.get_x(np_to, tx);
  interp.optclear_and_resize(npt2*npf2);
  for (Int i = 0; i < np_to; ++i) {
    Real fi[GLL::np_max];
    GLL::eval_lagrange_poly(np_from, fx, tx[i], fi);
    for (Int j = 0; j < np_to; ++j) {
      Real fj[GLL::np_max];
      GLL::eval_lagrange_poly(np_from, fx, tx[j], fj);
      for (Int ii = 0; ii < np_from; ++ii)
        for (Int ij = 0; ij < np_from; ++ij)
          interp[npf2*(np_to*i + j) + (np_from*ii + ij)] = fi[ii]*fj[ij];
    }
  }  
}

void IdemFv2Gll::build_npi_to_np_matrix (const Int np_from, const Int np_to,
                                         const Basis& b, Array<Real>& op) {
  build_interp_matrix(np_from, np_to, b, op);
}

void IdemFv2Gll::remapd (const Real* gll_metdet, const Real* fv_metdet,
                         const Real* p, Real* d) const {
  const Int nphys2 = nphys*nphys, npi2 = npi*npi, np2 = np*np;
  assert(npi <= np);
  assert(nphys2 <= npi2);
  Real wrk[GLL::np_max*GLL::np_max];
  for (Int i = 0; i < nphys2; ++i)
    wrk[i] = fv_metdet[i]*p[i];
  matvec(np2, nphys2, op_p_to_d.data(), wrk, d);
  for (Int i = 0; i < np2; ++i)
    d[i] /= gll_metdet[i];
}

void IdemFv2Gll
::reconstruct_nphys1 (const Real* metdet, Real* rho, Real* q) const {
  if (np <= np_nphys1) return;
  const Int np2 = np*np;
  for (Int i = 0; i < np2; ++i) rho[i] *= metdet[i];
  for (Int i = 0; i < np2; ++i) q[i] *= rho[i];
  reconstructd_nphys1(rho);
  reconstructd_nphys1(q);
  for (Int i = 0; i < np2; ++i) q[i] /= rho[i];
  for (Int i = 0; i < np2; ++i) rho[i] /= metdet[i];
}

void IdemFv2Gll::reconstructd_nphys1 (Real* rho) const {
  const Int npf = np_nphys1, npf2 = npf*npf, np2 = np*np;
  Real wrk[np_nphys1*np_nphys1];
  if (np_nphys1 == 3) {
    // Can't use this unless we write region limiter code. np_nphys1 2 and 3
    // have approx the same empirical OOA, ~1.9, but 3 gives more accuracy. So
    // we may want to impl this case if pg1 becomes important.
    assert(0);
    // Combine regions of np cell to give an npf cell such that a DSS isn't
    // needed.
    const auto get_range = [&] (const Int i) -> std::pair<Int,Int> {
      if (i == 0) return std::make_pair(0, 1);
      else if (i == 2) return std::make_pair(np-1, np);
      return std::make_pair(1, np-1);
    };
    for (Int i = 0; i < npf; ++i) {
      const auto ir = get_range(i);
      for (Int j = 0; j < npf; ++j) {
        const auto jr = get_range(j);
        Real accum = 0, den = 0;
        for (Int si = ir.first; si < ir.second; ++si)
          for (Int sj = jr.first; sj < jr.second; ++sj) {
            const Int idx = np*si + sj;
            accum += w_dd[idx]*rho[idx];
            den += w_dd[idx];
          }
        wrk[npf*i + j] = accum / den;
      }
    }
  } else {
    // Use just the corner points. The idea is that the non-corner points have
    // not interacted with the corner points. Thus, the corners behave as though
    // everything has been done on a GLL np=2 grid. In particular, this means
    // that using just the corner values now is still mass conserving.
    wrk[0] = rho[0];
    wrk[1] = rho[np-1];
    wrk[2] = rho[(np-1)*np];
    wrk[3] = rho[np2-1];
  }
  // Interpolate from npf to np. Because each edge is interpolated using only
  // its edge data, a DSS is not needed after this.
  for (Int i = 0; i < np2; ++i) {
    Real accum = 0;
    for (Int j = 0; j < npf2; ++j)
      accum += npi_to_np[npf2*i + j]*wrk[j];
    rho[i] = accum;
  }
}

void L2Fv2Gll::init (Int np, Int nphys, const Basis::ConstPtr& basis) {
  Remap::init(np, nphys, basis);
  npi = nphys == 1 ? 2 : nphys;
  init_matrices();
}

void L2Fv2Gll::init_matrices () {
  // LS problem for densities p and d is
  //     diag(gll_metdet) d = M_dd \ (M_ddi (M_didi \ (M_dip diag(fv_metdet) p)))
  // Form op_p_to_d = M_dd \ (M_ddi (M_didi \ M_dip)) so that
  //     diag(gll_metdet) d = op_p_to_d diag(fv_metdet) p.
  const Int nphys2 = nphys*nphys, npi2 = npi*npi, np2 = np*np;
  Array<Real> M_pdi(nphys2*npi2), M_didi(npi2*npi2), M_ddi(np2*npi2), M_dd(np2*np2),
    e(nphys2, 0), x(np2, 0), y(np2, 0);
  UniformNodeReduced b_npi;
  Basis::compute_integrals_over_subcells_2d(b_npi, npi, nphys, M_pdi.data());
  Basis::compute_mass_matrix_2d(b_npi, npi, b_npi, npi, M_didi.data());
  Basis::compute_mass_matrix_2d(*b, np, b_npi, npi, M_ddi.data());
  Basis::compute_mass_matrix_2d(*b, np, *b, np, M_dd.data());
  dpotrf('L', npi2, M_didi.data(), npi2);
  dpotrf('L', np2, M_dd.data(), np2);
  op_p_to_d.optclear_and_resize(np2*nphys2);
  for (Int i = 0; i < nphys2; ++i) {
    e[i] = 1;
    tmatvec(nphys2, npi2, M_pdi.data(), e.data(), x.data());
    dpotrs('L', npi2, 1, M_didi.data(), npi2, x.data(), npi2);
    matvec(np2, npi2, M_ddi.data(), x.data(), y.data());
    dpotrs('L', np2, 1, M_dd.data(), np2, y.data(), np2);
    for (Int j = 0; j < np2; ++j) op_p_to_d[nphys2*j + i] = y[j];
    e[i] = 0;
  }
}

void L2Fv2Gll::remapd (const Real* gll_metdet, const Real* fv_metdet,
                       const Real* p, Real* d) const {
  const Int nphys2 = nphys*nphys, npi2 = npi*npi, np2 = np*np;
  assert(npi <= np);
  assert(nphys2 <= npi2);
  Real wrk[GLL::np_max*GLL::np_max];
  for (Int i = 0; i < nphys2; ++i)
    wrk[i] = fv_metdet[i]*p[i];
  matvec(np2, nphys2, op_p_to_d.data(), wrk, d);
  for (Int i = 0; i < np2; ++i)
    d[i] /= gll_metdet[i];
}

void L2ExceptPerimFv2Gll::init (Int np, Int nphys, const Basis::ConstPtr& basis) {
  Remap::init(np, nphys, basis);
  npi = nphys == 1 ? 2 : nphys;
  init_matrices();
}

void L2ExceptPerimFv2Gll::init_matrices () {
  /* First remap problem, fv_metdet p = mp -> di. p is the set of perimeter
     subcell indices.
        min_di 1/2 di'M_didi di - di'M_dip mp
         st    M_dip(:,p)'di = M_pp(p,p) mp(p)  ; conserve mass in each perimeter subcell
               w_di'di = w_p'mp                 ; mass conservation over whole element
     Let Con = [M_dip(:,p)'; w_di'] and r = [mp(p); w_p'mp].
        Lag = 1/2 di'M_didi di - di'M_dip p + mu'(Con di - r)
        Lag_di = M_didi di - M_dip mp + Con mu = 0
        Lag_mu = Con'di - r = 0
     =>
        [M_didi Con] [di] = [M_dip mp] = [M_dip    ] [mp]
        [Con'   0  ] [mu]   [r       ]   [[p; w_p']]
     Second remap problem: di -> gll_metdet d = md
        min_md 1/2 md'M_dd md - md'M_ddi di
     => M_dd md = M_ddi di
   */
  const Int nphys2 = nphys*nphys, npi2 = npi*npi, np2 = np*np;
  UniformNodeReduced b_npi;
  Array<Real> CX(npi2*nphys2);
  {
    const bool fully_constrained = nphys == 2;
    const Int mass_constraint = fully_constrained ? 0 : 1;
    const Int ncon = 4*(nphys - 1) + mass_constraint;
    Array<Real> M_pdi(nphys2*npi2), M_didi(npi2*npi2), Con(npi2*ncon), D(ncon*nphys2, 0),
      rwrk(npi2*ncon);
    Array<Int> iwrk(ncon);
    Basis::compute_integrals_over_subcells_2d(b_npi, npi, nphys, M_pdi.data());
    Basis::compute_mass_matrix_2d(b_npi, npi, b_npi, npi, M_didi.data());
    // All but last rows of Con and D.
    Int k = 0;
    for (Int i = 0; i < nphys; ++i)
      for (Int j = 0; j < nphys; ++j) {
        // Want only perimeter subcells.
        if (i > 0 && i < nphys-1 && j > 0 && j < nphys-1) continue;
        // perimeter subcell index perim(k)
        const Int subcell = i*nphys + j;
        // M_pdi(perim(k),:)
        const Real* const M_pdi_row = &M_pdi[npi2*subcell];
        std::copy(M_pdi_row, M_pdi_row + npi2, &Con[npi2*k]);
        // perim index matrix
        D[ncon*subcell + k] = 4.0/nphys2;
        ++k;
      }
    assert(k == ncon - mass_constraint);
    if (mass_constraint) {
      // Final constraint matrix row.
      const Real* w;
      b_npi.get_w(npi, w);
      Real* const w_didi = &Con[npi2*(ncon-1)];
      for (Int i = 0; i < npi; ++i)
        for (Int j = 0; j < npi; ++j)
          w_didi[npi*i + j] = w[i]*w[j];
      // Final row of D.
      for (Int i = 0; i < nphys2; ++i)
        D[ncon*(i+1) - 1] = 4.0/nphys2;
    }
    std::copy(M_pdi.begin(), M_pdi.end(), CX.begin());
    const Int status =
      solve_kkt(npi2, ncon, nphys2, M_didi.data(), Con.data(), CX.data(), D.data(),
                rwrk.data(), iwrk.data());
    SIQK_THROW_IF(status != 0, "KKT system in L2ExceptPerimFv2Gll could not be solved.");
  }
  {
    Array<Real> M_ddi(np2*npi2), M_dd(np2*np2), x(np2*nphys2, 0);
    Basis::compute_mass_matrix_2d(*b, np, b_npi, npi, M_ddi.data());
    Basis::compute_mass_matrix_2d(*b, np, *b, np, M_dd.data());
    matmult_rcc(np2, nphys2, npi2, 1, M_ddi.data(), CX.data(), 0, x.data());
    dpotrf('L', np2, M_dd.data(), np2);
    dpotrs('L', np2, nphys2, M_dd.data(), np2, x.data(), np2);
    op_p_to_d.optclear_and_resize(np2*nphys2);
    for (Int j = 0; j < np2; ++j)
      for (Int i = 0; i < nphys2; ++i)
        op_p_to_d[nphys2*j + i] = x[np2*i + j];
  }
}

void L2ExceptPerimFv2Gll
::remapd (const Real* gll_metdet, const Real* fv_metdet, const Real* p, Real* d) const {
  const Int nphys2 = nphys*nphys, npi2 = npi*npi, np2 = np*np;
  assert(npi <= np);
  assert(nphys2 <= npi2);
  Real wrk[GLL::np_max*GLL::np_max];
  for (Int i = 0; i < nphys2; ++i)
    wrk[i] = fv_metdet[i]*p[i];
  matvec(np2, nphys2, op_p_to_d.data(), wrk, d);
  for (Int i = 0; i < np2; ++i)
    d[i] /= gll_metdet[i];
}

void ElemLclReconFv2Gll::init (Int np, Int nphys, const Basis::ConstPtr& basis) {
  Remap::init(np, nphys, basis);
  init_matrices();
}

static bool
compute_integrals_over_subcells_2d (const Basis& b, const Int npi, const Int npj,
                                    const Int nfi, const Int nfj, Real* M,
                                    const bool test) {
  const auto subcell_idx = [&] (const Int nf, const Real r) -> Int {
    const Real p = (1 + r)/2;
    return std::max(0, std::min(nf-1, Int(nf*p)));
  };

  // The common refinement here and later is needed only if the Basis b is not
  // smooth, e.g. an Islet basis, but is fine to use (of course) for smooth
  // bases.
  std::vector<Real> xcom, ycom;
  if ( ! Basis::calc_common_refinement(b, npi, nfi, ycom, test) ||
       ! Basis::calc_common_refinement(b, npj, nfj, xcom, test)) {
    assert(0);
    return false; 
  }
  const Int nxcom = xcom.size(), nycom = ycom.size(), npij = npi*npj, nfij = nfi*nfj;

  // Quadrature data.
  Int qn = (b.max_degree(std::max(npi, npj)) + 4)/2;
  if (qn > 13 && qn < 16) qn = 16;
  const Real* qx, * qwt;
  GLL q_gll; // GLL here rather than the basis b/c this is for the quadrature rule
  q_gll.get_coef(qn, qx, qwt);

  for (Int i = 0; i < npij*nfij; ++i) M[i] = 0;

  Real vx[GLL::np_max], vy[GLL::np_max];
  for (Int iy = 0; iy < nycom-1; ++iy) {
    const Real yb = ycom[iy], ye = ycom[iy+1];
    const Int ify = subcell_idx(nfi, 0.5*(yb + ye));
    for (Int ix = 0; ix < nxcom-1; ++ix) {
      const Real xb = xcom[ix], xe = xcom[ix+1];
      const Int ifx = subcell_idx(nfj, 0.5*(xb + xe));
      Real* Msub = M + npij*(ify*nfj + ifx);
      const Real fac = 0.25*(ye - yb)*(xe - xb);
      for (Int jy = 0; jy < qn; ++jy) {
        Real alpha = 0.5*(qx[jy] + 1);
        const Real y = (1 - alpha)*yb + alpha*ye;
        b.eval(npi, y, vy);
        for (Int jx = 0; jx < qn; ++jx) {
          alpha = 0.5*(qx[jx] + 1);
          const Real x = (1 - alpha)*xb + alpha*xe;
          b.eval(npj, x, vx);
          const Real fac1 = fac*qwt[jy]*qwt[jx];
          for (Int i = 0; i < npi; ++i)
            for (Int j = 0; j < npj; ++j)
              Msub[npj*i + j] += fac1*vy[i]*vx[j];
        }
      }
    }
  }

  return true;
}

struct MixedOpData {
  Int ni, nj;
  Array<Real> R, tau, wrk;
};

static bool make_mixed_op (const Int ni, const Int nj, const Basis& b, MixedOpData& op) {
  const Int nij = ni*nj, nij2 = nij*nij;
  op.ni = ni; op.nj = nj;
  op.R.optclear_and_resize(nij2);
  op.wrk.optclear_and_resize(nij2);
  bool ok = compute_integrals_over_subcells_2d(b, ni, nj, ni, nj, op.R.data(), true);
  if ( ! ok) return false;
  op.tau.optclear_and_resize(nij2);
  Array<Real> wrk(nij2);
  ok = dgeqrf(nij, nij, op.R.data(), nij, op.tau.data(), op.wrk.data(), nij2) == 0;
  return ok;
}

static bool apply_mixed_op (const MixedOpData& op, const Int nrhs, Real* x) {
  const Int n = op.ni*op.nj;
  for (Int i = 0; i < n; ++i) x[i] *= 4.0/n;
  const bool ok =
    (dtrsm('l', 'u', 't', 'n', n, nrhs, 1.0, op.R.data(), n, x, n) == 0 &&
     dormqr('l', 'n', n, nrhs, n, op.R.data(), n, op.tau.data(), x, n,
            op.wrk.data(), n) == 0);
  return ok;
}

struct Panel {
  static constexpr Int edge_np = 3, interior_np = 3,
    cap = interior_np*interior_np;
  Int ni, nj;
  Real v[cap];

  // 1D coordinates.
  static void set_coord (const Int nf, const Int j,      // subcell index in element
                         const Real& xb, const Real& xe, // interval w.r.t. element
                         Real& xb_pnl, Real& xe_pnl) {   // interval w.r.t. panel
    assert(edge_np == 2 || edge_np == 3);
    
    Real p[2]; // full extent of subcell w.r.t. panel
    if (j == 0 || j == nf-1) {
      if (edge_np == 2) {
        if (j == 0) { p[0] = -1; p[1] = 0; }
        else { p[0] = 0; p[1] = 1; }
      } else {
        static const Real panel_middle = 1.0/edge_np;
        if (j == 0) { p[0] = -1; p[1] = -panel_middle; }
        else { p[0] = panel_middle; p[1] = 1; }
      }
    } else {
      static const Real panel_middle = 1.0/interior_np;
      p[0] = -panel_middle; p[1] = panel_middle;
    }
    Real e[2]; // full extent of subcell w.r.t. element
    e[0] = 2*Real(j)/nf - 1; e[1] = 2*Real(j+1)/nf - 1;
    assert(xb >= e[0] && xe <= e[1]);
    Real alpha = (xb - e[0])/(e[1] - e[0]);
    xb_pnl = (1 - alpha)*p[0] + alpha*p[1];
    alpha = (xe - e[0])/(e[1] - e[0]);
    xe_pnl = (1 - alpha)*p[0] + alpha*p[1];
  }

  static Int subcell_idx (const Int nf, const Real r) {
    const Real p = (1 + r)/2;
    return std::max(0, std::min(nf-1, Int(nf*p)));
  }
};

// Apply the IdemFv2Gll method to each panel to reconstruct, for each panel, one
// subcell. Then L2 project this reconstruction to the np basis.
void ElemLclReconFv2Gll::init_matrices () {
  const Int nf = nphys, nf2 = square(nf), np2 = square(np);
  bool ok;

  // 1D pg-b common refinement
  std::vector<Real> xcom;
  ok = Basis::calc_common_refinement(*b, np, nf, xcom);
  assert(ok);
  const Int ncom = xcom.size();

  GLL p_gll; // GLL and not basis because this is the panel basis
  MixedOpData op_corner, op_edge_jint, op_edge_iint, op_int;
  { // Panel operators; fast index is second argument.
    const Int enp = Panel::edge_np, inp = Panel::interior_np;
    ok = (make_mixed_op(enp, enp, p_gll, op_corner) &&
          make_mixed_op(enp, inp, p_gll, op_edge_jint) &&
          make_mixed_op(inp, enp, p_gll, op_edge_iint) &&
          make_mixed_op(inp, inp, p_gll, op_int));
    assert(ok);
  }

  const Int edge_np = Panel::edge_np, interior_np = Panel::interior_np;
  Int qn = (p_gll.max_degree(std::max(edge_np, interior_np)) +
            b->max_degree(np) + 4)/2;
  if (qn > 13 && qn < 16) qn = 16;
  const Real* qx, * qwt;
  GLL q_gll; // GLL and not basis b/c this is the quadrature rule
  q_gll.get_coef(qn, qx, qwt);

  Array<Real> M_mix(np2*nf2);
  { // Mixed mass matrix; fast index is np-basis.
    Array<Real> ei(nf2, 0);
    std::vector<Panel> panels(nf2);
    for (Int dof_pg = 0; dof_pg < nf2; ++dof_pg) {
      ei[dof_pg] = 1;
      // Form a reconstruction of the ei function in each subcell, represented
      // as panels.
      for (Int sci = 0; sci < nf; ++sci) {
        const bool iedge = sci == 0 || sci == nf-1;
        const Int ihalo = iedge ? edge_np-1 : interior_np/2;
        for (Int scj = 0; scj < nf; ++scj) {
          const bool jedge = scj == 0 || scj == nf-1;
          const Int jhalo = jedge ? edge_np-1 : interior_np/2;
          auto& pnl = panels[sci*nf + scj];
          pnl.ni = iedge ? Panel::edge_np : Panel::interior_np;
          pnl.nj = jedge ? Panel::edge_np : Panel::interior_np;
          const Int ndof_pg = pnl.ni*pnl.nj;
          // Get relevant p entries.
          Int k = 0;
          for (Int i = 0; i < nf; ++i)
            for (Int j = 0; j < nf; ++j) {
              if (std::abs(i - sci) > ihalo || std::abs(j - scj) > jhalo) continue;
              pnl.v[k] = ei[nf*i + j];
              ++k;
            }
          assert(k == ndof_pg);
          // Solve for d entries.
          const bool edge = iedge || jedge, corner = iedge && jedge;
          const auto& idem_op = ( ! edge ? op_int :
                                  (corner ? op_corner :
                                   (jedge ? op_edge_iint : op_edge_jint)));
          assert(square(ndof_pg) == (Int)idem_op.R.size());
          ok = apply_mixed_op(idem_op, 1, pnl.v);
          assert(ok);
        }
      }
      ei[dof_pg] = 0;
      // Compute entries of mixed mass matrix using basis and panels.
      Real vx[GLL::np_max], vy[GLL::np_max],
        vx_pnl[Panel::interior_np], vy_pnl[Panel::interior_np];
      Real* const M_mix_col = &M_mix[np2*dof_pg];
      for (Int i = 0; i < np2; ++i) M_mix_col[i] = 0;
      for (Int iy = 0; iy < ncom-1; ++iy) { // iterate over common refinement
        const Real yb = xcom[iy], ye = xcom[iy+1]; // coordinate in element
        const Int sci = Panel::subcell_idx(nf, 0.5*(yb + ye)); // subcell index
        Real yb_pnl, ye_pnl; // coord w.r.t. panel
        Panel::set_coord(nf, sci, yb, ye, yb_pnl, ye_pnl);
        for (Int ix = 0; ix < ncom-1; ++ix) {
          const Real xb = xcom[ix], xe = xcom[ix+1];
          const Int scj = Panel::subcell_idx(nf, 0.5*(xb + xe));
          Real xb_pnl, xe_pnl;
          Panel::set_coord(nf, scj, xb, xe, xb_pnl, xe_pnl);
          // Integrate over this facet.
          const Real fac = 0.25*(ye - yb)*(xe - xb);
          const auto& pnl = panels[sci*nf + scj];
          for (Int jy = 0; jy < qn; ++jy) {
            Real alpha = 0.5*(qx[jy] + 1);
            const Real y_pnl = (1 - alpha)*yb_pnl + alpha*ye_pnl;
            p_gll.eval(pnl.ni, y_pnl, vy_pnl);
            const Real y = (1 - alpha)*yb + alpha*ye;
            b->eval(np, y, vy);
            for (Int jx = 0; jx < qn; ++jx) {
              alpha = 0.5*(qx[jx] + 1);
              const Real x_pnl = (1 - alpha)*xb_pnl + alpha*xe_pnl;
              p_gll.eval(pnl.nj, x_pnl, vx_pnl);
              const Real x = (1 - alpha)*xb + alpha*xe;
              b->eval(np, x, vx);
              // Panel reconstruction value.
              Real f_pnl = 0;
              for (Int i = 0; i < pnl.ni; ++i)
                for (Int j = 0; j < pnl.nj; ++j)
                  f_pnl += vy_pnl[i]*vx_pnl[j]*pnl.v[pnl.nj*i + j];
              // Integrate against basis functions.
              const Real fac1 = fac*qwt[jy]*qwt[jx];
              for (Int i = 0; i < np; ++i)
                for (Int j = 0; j < np; ++j)
                  M_mix_col[np*i + j] += fac1*f_pnl*vy[i]*vx[j];
            }
          }
        }
      }
    }
  }
  
  // Form rest of L2 op.
  Array<Real> M_dd(np2*np2);
  Basis::compute_mass_matrix_2d(*b, np, *b, np, M_dd.data());
  dpotrf('L', np2, M_dd.data(), np2);
  dpotrs('L', np2, nf2, M_dd.data(), np2, M_mix.data(), np2);
  op_p_to_d.optclear_and_resize(np2*nf2);
  for (Int i = 0; i < np2; ++i)
    for (Int j = 0; j < nf2; ++j)
      op_p_to_d[nf2*i + j] = M_mix[np2*j + i];
}

void ElemLclReconFv2Gll
::remapd (const Real* gll_metdet, const Real* fv_metdet, const Real* p, Real* d) const {
  const Int nphys2 = nphys*nphys, np2 = np*np;
  Real wrk[GLL::np_max*GLL::np_max];
  for (Int i = 0; i < nphys2; ++i)
    wrk[i] = fv_metdet[i]*p[i];
  matvec(np2, nphys2, op_p_to_d.data(), wrk, d);
  for (Int i = 0; i < np2; ++i)
    d[i] /= gll_metdet[i];
}

void MeshData::init (const Mesh& m, const Real* jacobian_gll) {
  nelem = nslices(m.geo_c2n);
  ncgll = nslices(m.cgll_p);
  ndgll = nslices(m.dglln2cglln);
  if (jacobian_gll) {
    resize(gll_metdet, ndgll);
    std::copy(jacobian_gll, jacobian_gll + ndgll, gll_metdet.data());
  } else {
    calc_node_jacobians(m, m.geo_p, gll_metdet);
  }
  calc_gll_basis_function_integrals(m, *m.basis, spheremp);
  d2cer = std::make_shared<D2Cer>(m.dglln2cglln, spheremp);
  wrk.optclear_and_resize(ndgll);
  mesh::get_adjacent_cells(m.geo_c2n, geo_c2cnbrs_ptr, geo_c2cnbrs);
}

void Remap::print () const {
  printf("w_dd (%d):", np);
  for (Int i = 0; i < np*np; ++i) printf("%12.4e", w_dd[i]);
  printf("\n");
}

void Gll2Fv::print () const {
  printf("Gll2Fv:\n");
  Remap::print();
  printf("M_pp (%d):", nphys);
  for (Int i = 0; i < nphys*nphys; ++i) printf("%12.4e", M_pp[i]);
  printf("\nM_dp (%d,%d):", np, nphys);
  for (Int i = 0; i < np*np*nphys*nphys; ++i) printf("%12.4e", M_dp[i]);
  printf("\n");
}

// We have metdet for the GLL mesh. Compute it for the FV mesh by remapping the
// unit density field.
void init_metdet (MeshData& md, const Gll2Fv& pg) {
  const Int nphys = pg.get_nphys(), nphys2 = nphys*nphys,
    np = pg.get_np(), np2 = np*np;
  md.fv_metdet.optclear_and_resize(md.nelem*nphys2);
  Array<Real> ones(std::max(np2, nphys2), 1);
  for (Int ie = 0; ie < md.nelem; ++ie)
    pg.remapd(&md.gll_metdet[ie*np2],
              ones.data(), // pretend metdet on the FV mesh is 1
              ones.data(), // unit density field
              // compute metdet on the FV mesh
              &md.fv_metdet[ie*nphys2]);
}

PhysgridOps::PhysgridOps (const Mesh& mesh, const Int np_, const Int nphys_,
                          const Fv2Gll::Type::Enum fv2gll_type,
                          const Real* jacobian_gll)
  : np(np_), nphys(nphys_)
{
  SIQK_THROW_IF(nphys > np, "Physgrid remap supports only nphys <= np.");
  mesh_data = std::make_shared<MeshData>(mesh, jacobian_gll);
  basis = mesh.basis;
  gll2fv = std::make_shared<Gll2Fv>(np, nphys, basis);
  fv2gll = Fv2Gll::create(fv2gll_type, np, nphys, basis);
  init_metdet(*mesh_data, *gll2fv);
}

vis::VisWriter::Ptr make_pg_vis (const std::string& fname_prefix, Int ne, Int nphys,
                                 Int res, bool nonuniform) {
  AVec3s pg_p;
  AIdxs pg_c2n;
  const auto basis = Basis::create(Basis::Type::uniform_offset_nodal);
  mesh::make_cubedsphere_subcell_mesh(ne, nphys+1, *basis, pg_p, pg_c2n);
  if (nonuniform) mesh::make_nonuniform(pg_p);
  const auto op_pg = std::make_shared<vis::PhysgridToLatLon>(pg_p, pg_c2n, 2*res+1, 4*res+1);
  return std::make_shared<vis::VisWriter>(op_pg, fname_prefix + "_pg.bin");
}

} // namespace pg
