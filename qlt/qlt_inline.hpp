#ifndef INCLUDE_QLT_INLINE_HPP
#define INCLUDE_QLT_INLINE_HPP

#include <cassert>

namespace qlt {

template <typename ES> KOKKOS_INLINE_FUNCTION
void QLT<ES>::set_rho (const Int& lclcellidx, const Real& rhom) {
  const Int ndps = md_.a_d.prob2bl2r[md_.nprobtypes];
  bd_.l2r_data(ndps*lclcellidx) = rhom;  
}

template <typename ES> KOKKOS_INLINE_FUNCTION
void QLT<ES>::set_Q (const Int& lclcellidx, const Int& tracer_idx,
                     const Real& Qm,
                     const Real& Qm_min, const Real& Qm_max,
                     const Real Qm_prev) {
  const Int ndps = md_.a_d.prob2bl2r[md_.nprobtypes];
  Real* bd; {
    const Int bdi = md_.a_d.trcr2bl2r(tracer_idx);
    bd = &bd_.l2r_data(ndps*lclcellidx + bdi);
  }
  bd[1] = Qm;
  {
    const Int problem_type = md_.a_d.trcr2prob(tracer_idx);
    if (problem_type & ProblemType::shapepreserve) {
      bd[0] = Qm_min;
      bd[2] = Qm_max;
    } else if (problem_type & ProblemType::consistent) {
      const Real rhom = bd_.l2r_data(ndps*lclcellidx);
      bd[0] = Qm_min / rhom;
      bd[2] = Qm_max / rhom;
    } else {
      Kokkos::abort("set_Q: invalid problem_type.");
    }
    if (problem_type & ProblemType::conserve) {
      if (Qm_prev < 0) Kokkos::abort("Qm_prev was not provided to set_Q.");
      bd[3] = Qm_prev;
    }
  }
}

template <typename ES> KOKKOS_INLINE_FUNCTION
Real QLT<ES>::get_Q (const Int& lclcellidx, const Int& tracer_idx) {
  const Int ndps = md_.a_d.prob2br2l[md_.nprobtypes];
  const Int bdi = md_.a_d.trcr2br2l(tracer_idx);
  return bd_.r2l_data(ndps*lclcellidx + bdi);
}

namespace impl {
// GPU-friendly replacements for std::min/max.
template <typename T> KOKKOS_INLINE_FUNCTION
const T& min (const T& a, const T& b) { return a < b ? a : b; }
template <typename T> KOKKOS_INLINE_FUNCTION
const T& max (const T& a, const T& b) { return a > b ? a : b; }
}

namespace slv {
KOKKOS_INLINE_FUNCTION
Real get_xbd (const Real* xbd, const Int i, const bool xbds_scalar)
{ return xbds_scalar ? *xbd : xbd[i]; }

KOKKOS_INLINE_FUNCTION
bool is_inside (const Real xi, const Real* xlo, const Real* xhi, const Int i,
                const bool xbds_scalar) {
  return (xi > get_xbd(xlo, i, xbds_scalar) &&
          xi < get_xbd(xhi, i, xbds_scalar));
}

KOKKOS_INLINE_FUNCTION
bool is_outside (const Real xi, const Real* xlo, const Real* xhi, const Int i,
                 const bool xbds_scalar) {
  return (xi < get_xbd(xlo, i, xbds_scalar) ||
          xi > get_xbd(xhi, i, xbds_scalar));
}

KOKKOS_INLINE_FUNCTION
Real calc_r_tol (const Real b, const Real* a, const Real* y, const Int n) {
  Real ab = std::abs(b);
  for (Int i = 0; i < n; ++i) ab = std::max(ab, std::abs(a[i]*y[i]));
  return 1e1*std::numeric_limits<Real>::epsilon()*std::abs(ab);
}

KOKKOS_INLINE_FUNCTION
void calc_r (const Int n, const Real* w, const Real* a, const Real b,
             const Real* xlo, const Real* xhi, const bool xbds_scalar,
             const Real* y, const Real& lambda, Real* x, Real& r, Real& r_lambda) {
  r = 0;
  r_lambda = 0;
  for (Int i = 0; i < n; ++i) {
    const Real q = a[i]/w[i];
    const Real x_trial = y[i] + lambda*q;
    Real xtmp;
    if (x_trial < (xtmp = get_xbd(xlo, i, xbds_scalar)))
      x[i] = xtmp;
    else if (x_trial > (xtmp = get_xbd(xhi, i, xbds_scalar)))
      x[i] = xtmp;
    else {
      x[i] = x_trial;
      r_lambda += a[i]*q;
    }
    r += a[i]*x[i];
  }
  r -= b;
}

// Solve
//     min_x sum_i w(i) (x(i) - y(i))^2
//      st   a' x = b
//           xlo <= x <= xhi.
// This function assumes w > 0 to save a few operations. Return 0 on success and
// x == y, 1 on success and x != y, -1 if infeasible, -2 if max_its hit with no
// solution. See Section 3 of Bochev, Ridzal, Shashkov, Fast optimization-based
// conservative remap of scalar fields through aggregate mass transfer. lambda
// is used in check_1eq_bc_qp_foc.
//todo 2D version of this function that takes advantage of 2D.
KOKKOS_INLINE_FUNCTION
Int solve_1eq_bc_qp (const Int n, const Real* w, const Real* a, const Real b,
                     const Real* xlo, const Real* xhi, const bool xbds_scalar,
                     const Real* y, Real* x, const Int max_its = 100) {
  const Real r_tol = calc_r_tol(b, a, y, n);

  { // Check for a quick exit.
    bool all_in = true;
    Real r = 0;
    for (Int i = 0; i < n; ++i) {
      if (is_outside(x[i], xlo, xhi, i, xbds_scalar)) {
        all_in = false;
        break;
      }
      r += a[i]*x[i];
    }
    if (all_in) {
      r -= b;
      if (std::abs(r) <= r_tol)
        return 0;
    }
  }

  { // Eval r at end points to check for feasibility, and also possibly a quick
    // exit on a common case.
    Real r = -b;
    for (Int i = 0; i < n; ++i) {
      x[i] = get_xbd(xlo, i, xbds_scalar);
      r += a[i]*x[i];
    }
    if (std::abs(r) <= r_tol) return 1;
    if (r > 0) return -1;
    r = -b;
    for (Int i = 0; i < n; ++i) {
      x[i] = get_xbd(xhi, i, xbds_scalar);
      r += a[i]*x[i];
    }
    if (std::abs(r) <= r_tol) return 1;
    if (r < 0) return -1;
  }

  { // Check for a quick exit: the bounds are so tight that the midpoint of the
    // box satisfies r_tol.
    Real r = -b;
    for (Int i = 0; i < n; ++i) {
      x[i] = 0.5*(get_xbd(xlo, i, xbds_scalar) + get_xbd(xhi, i, xbds_scalar));
      r += a[i]*x[i];
    }
    if (std::abs(r) <= r_tol) return 1;
  }

  const Real wall_dist = 1e-3;

  // Get lambda endpoints.
  Real lamlo = 0, lamhi = 0;
  for (Int i = 0; i < n; ++i) {
    const Real rq = w[i]/a[i];
    const Real lamlo_i = rq*(get_xbd(xlo, i, xbds_scalar) - y[i]);
    const Real lamhi_i = rq*(get_xbd(xhi, i, xbds_scalar) - y[i]);
    if (i == 0) {
      lamlo = lamlo_i;
      lamhi = lamhi_i;
    } else {
      lamlo = impl::min(lamlo, lamlo_i);
      lamhi = impl::max(lamhi, lamhi_i);
    }
  }
  const Real lamlo_feas = lamlo, lamhi_feas = lamhi;
  Real lambda = lamlo <= 0 && lamhi >= 0 ? 0 : lamlo;

  Int info = -2;

  // Bisection-safeguarded Newton iteration for r(lambda) = 0.
  bool prev_step_bisect = false;
  Int nbisect = 0;
  for (Int iteration = 0; iteration < max_its; ++iteration) {
    // Compute x, r, r_lambda.
    Real r, r_lambda;
    calc_r(n, w, a, b, xlo, xhi, xbds_scalar, y, lambda, x, r, r_lambda);
    // Is r(lambda) - b sufficiently == 0?
    if (std::abs(r) <= r_tol) {
      info = 1;
      break;
    }
    // Check if the lambda bounds are too close.
    if (nbisect > 64) {
      if (lamhi == lamhi_feas || lamlo == lamlo_feas) {
        // r isn't small enough and one lambda bound is on the feasibility
        // limit. The QP must not be feasible.
        info = -1;
        break;
      }
      info = 1;
      break;
    }
    // Adjust lambda bounds.
    if (r > 0)
      lamhi = lambda;
    else
      lamlo = lambda;
    if (r_lambda != 0) {
      // Newton step.
      lambda -= r/r_lambda;
    } else {
      // Force bisection.
      lambda = lamlo;
    }
    // Safeguard. The wall distance check assures progress, but use it only
    // every other potential bisection.
    const Real D = prev_step_bisect ? 0 : wall_dist*(lamhi - lamlo);
    if (lambda - lamlo < D || lamhi - lambda < D) {
      lambda = 0.5*(lamlo + lamhi);
      ++nbisect;
      prev_step_bisect = true;
    } else {
      prev_step_bisect = false;
    }
  }

  return info;
}

KOKKOS_INLINE_FUNCTION
void r2l_nl_adjust_bounds (Real Qm_bnd[2], const Real rhom[2], Real Qm_extra) {
  Real q[2];
  for (Int i = 0; i < 2; ++i) q[i] = Qm_bnd[i] / rhom[i];
  if (Qm_extra < 0) {
    Int i0, i1;
    if (q[0] >= q[1]) { i0 = 0; i1 = 1; } else { i0 = 1; i1 = 0; }
    const Real Qm_gap = (q[i1] - q[i0])*rhom[i0];
    if (Qm_gap <= Qm_extra) {
      Qm_bnd[i0] += Qm_extra;
      return;
    }
  } else {
    Int i0, i1;
    if (q[0] <= q[1]) { i0 = 0; i1 = 1; } else { i0 = 1; i1 = 0; }
    const Real Qm_gap = (q[i1] - q[i0])*rhom[i0];
    if (Qm_gap >= Qm_extra) {
      Qm_bnd[i0] += Qm_extra;
      return;
    }
  }
  { // Have to adjust both. Adjust so that the q bounds are the same. This
    // procedure assures that as long as rhom is conservative, then the
    // adjustment never pushes q_{min,max} out of the safety bounds.
    const Real Qm_tot = Qm_bnd[0] + Qm_bnd[1] + Qm_extra;
    const Real rhom_tot = rhom[0] + rhom[1];
    const Real q_tot = Qm_tot / rhom_tot;
    for (Int i = 0; i < 2; ++i)
      Qm_bnd[i] = q_tot*rhom[i];
  }
}

KOKKOS_INLINE_FUNCTION
void r2l_l_adjust_bounds (const Int np, Real* q_min, Real* q_max, const Real* rhom,
                          Real Qm_extra) {
  assert(0); // Not used right now, but want to eventually. Need to do some more analysis.
  static constexpr int max_np = 16;
  Real* const q_bnd = Qm_extra < 0 ? q_min : q_max;
  // Try solving a QP that adjusts a q bound.
  Real Qm = Qm_extra;
  Real w[max_np], q_bnd_min[max_np], q_bnd_max[max_np], q_bnd_orig[max_np];
  q_bnd_min[0] = q_min[0];
  q_bnd_max[0] = q_max[0];
  for (Int i = 0; i < np; ++i) {
    const Real rhomi = rhom[i];
    Qm += q_bnd[i]*rhomi;
    q_bnd_orig[i] = q_bnd[i];
    w[i] = rhomi;
    if (Qm_extra < 0) {
      q_bnd_min[0] = impl::min(q_bnd_min[0], q_min[i]);
      q_bnd_max[i] = q_max[i];
    } else {
      q_bnd_min[i] = q_min[i];
      q_bnd_max[0] = impl::max(q_bnd_max[0], q_max[i]);
    }
  }
  if (Qm_extra < 0)
    for (Int i = 1; i < np; ++i) q_bnd_min[i] = q_bnd_min[0];
  else
    for (Int i = 1; i < np; ++i) q_bnd_max[i] = q_bnd_max[0];
  // Check for feasibility.
  bool feasible; {
    Real Qm_lo = 0, Qm_hi = 0;
    for (Int i = 0; i < np; ++i) {
      Qm_lo += q_bnd_min[i]*w[i];
      Qm_hi += q_bnd_max[i]*w[i];
    }
    feasible = Qm_lo <= Qm && Qm <= Qm_hi;
  }
  if (feasible) {
    solve_1eq_bc_qp(np, w, w, Qm, q_bnd_min, q_bnd_max, false, q_bnd_orig, q_bnd);
  } else {
    // The QP isn't feasible, so set the bound to a constant.
    Real rhom_tot = 0, Qm_tot = Qm_extra;
    for (Int i = 0; i < np; ++i) {
      const Real rhomi = rhom[i];
      rhom_tot += rhomi;
      Qm_tot += q_bnd_orig[i]*rhomi;
    }
    const Real q_tot = Qm_tot / rhom_tot;
    for (Int i = 0; i < np; ++i)
      q_bnd[i] = q_tot;    
    //return;
    // Assert that this constant is outside of all previous bound values. That's
    // why the QP wasn't feasible.
    if (Qm_extra < 0)
      for (Int i = 0; i < np; ++i)
        assert(q_tot <= q_bnd_orig[i]);
    else
      for (Int i = 0; i < np; ++i)
        assert(q_tot >= q_bnd_orig[i]);
  }
}

KOKKOS_INLINE_FUNCTION
void solve_node_problem (const Real& rhom, const Real* pd, const Real& Qm,
                         const Real& rhom0, const Real* k0d, Real& Qm0,
                         const Real& rhom1, const Real* k1d, Real& Qm1) {
  Real Qm_min_kids [] = {k0d[0], k1d[0]};
  Real Qm_orig_kids[] = {k0d[1], k1d[1]};
  Real Qm_max_kids [] = {k0d[2], k1d[2]};
  { // Set the target values so that mass gets redistributed in a relative sense
    // rather than absolute. If a kid doesn't have much mass, don't give it too
    // much.
    const Real Qm_orig = pd[1], Qm_extra = Qm - Qm_orig;
    if (Qm_orig != 0)
      for (Int i = 0; i < 2; ++i)
        Qm_orig_kids[i] += (Qm_orig_kids[i] / Qm_orig) * Qm_extra;
  }
  { // The ideal problem is not assuredly feasible. Test for feasibility. If not
    // feasible, adjust bounds to solve the safety problem, which is assuredly
    // feasible if the total density field rho is mass conserving (Q doesn't
    // have to be mass conserving, of course; achieving mass conservation is one
    // use for QLT).
    const Real Qm_min = pd[0], Qm_max = pd[2];
    const bool lo = Qm < Qm_min, hi = Qm > Qm_max;
    if (lo || hi) {
      const Real rhom_kids[] = {rhom0, rhom1};
      r2l_nl_adjust_bounds(lo ? Qm_min_kids : Qm_max_kids,
                           rhom_kids,
                           Qm - (lo ? Qm_min : Qm_max));
    }
  }
  { // Solve the node's QP.
    static const Real ones[] = {1, 1};
    Real Qm_kids[2] = {k0d[1], k1d[1]};
    solve_1eq_bc_qp(2, ones, ones, Qm, Qm_min_kids, Qm_max_kids, false, Qm_orig_kids,
                    Qm_kids);
    Qm0 = Qm_kids[0];
    Qm1 = Qm_kids[1];
  }
}
} // namespace slv

template <typename ES> KOKKOS_INLINE_FUNCTION
void QLT<ES>::solve_node_problem (const Int problem_type,
                                  const Real& rhom, const Real* pd, const Real& Qm,
                                  const Real& rhom0, const Real* k0d, Real& Qm0,
                                  const Real& rhom1, const Real* k1d, Real& Qm1) {
  if ( ! (problem_type & ProblemType::shapepreserve)) {      
    Real mpd[3], mk0d[3], mk1d[3];
    mpd[0]  = pd [0]*rhom ; mpd [1] = pd[1] ; mpd [2] = pd [2]*rhom ;
    mk0d[0] = k0d[0]*rhom0; mk0d[1] = k0d[1]; mk0d[2] = k0d[2]*rhom0;
    mk1d[0] = k1d[0]*rhom1; mk1d[1] = k1d[1]; mk1d[2] = k1d[2]*rhom1;
    slv::solve_node_problem(rhom, mpd, Qm, rhom0, mk0d, Qm0, rhom1, mk1d, Qm1);
    return;
  }
  slv::solve_node_problem(rhom, pd, Qm, rhom0, k0d, Qm0, rhom1, k1d, Qm1);
}

} // namespace qlt

#endif
