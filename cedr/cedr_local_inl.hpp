#ifndef INCLUDE_CEDR_LOCAL_INL_HPP
#define INCLUDE_CEDR_LOCAL_INL_HPP

#include "cedr_util.hpp"

namespace cedr {
namespace local {

KOKKOS_INLINE_FUNCTION
bool is_inside (const Real xi, const Real* xlo, const Real* xhi, const Int i) {
  return (xi > xlo[i] && xi < xhi[i]);
}

KOKKOS_INLINE_FUNCTION
bool is_outside (const Real xi, const Real* xlo, const Real* xhi, const Int i) {
  return (xi < xlo[i] || xi > xhi[i]);
}

template <typename T>
KOKKOS_INLINE_FUNCTION
void sort4 (T* x) {
  T buf[4];
  for (Int i = 0; i < 2; ++i) {
    const Int j = 2*i;
    if (x[j] <= x[j+1]) { buf[j] = x[j]; buf[j+1] = x[j+1]; }
    else                { buf[j] = x[j+1]; buf[j+1] = x[j]; }
  }
  Int p0 = 0, p1 = 2;
  for (Int i = 0; i < 4; ++i)
    x[i] = (p1 >= 4 || (p0 < 2 && buf[p0] <= buf[p1]) ?
            buf[p0++] :
            buf[p1++]);
  cedr_assert(p0 == 2 && p1 == 4);
}

// 2D special case for efficiency.
KOKKOS_INLINE_FUNCTION
Int solve_1eq_bc_qp_2d (const Int n, const Real* w, const Real* a, const Real b,
                        const Real* xlo, const Real* xhi, 
                        const Real* y, Real* x) {
  cedr_assert(n == 2);

  { // Check if the optimal point ignoring bound constraints is in bounds.
    Real q[2], qsum = 0;
    for (int i = 0; i < 2; ++i) {
      q[i] = a[i]/w[i];
      qsum += q[i];
    }
    Real dm = b;
    for (int i = 0; i < 2; ++i)
      dm -= a[i]*y[i];
    bool good = true;
    for (int i = 0; i < 2; ++i) {
      x[i] = y[i] + dm*(q[i]/qsum);
      if (is_outside(x[i], xlo, xhi, i)) {
        good = false;
        break;
      }
    }
    if (good) return dm == 0 ? 0 : 1;
  }

  // Solve for intersection of a'x = b, given by the parameterized line
  //     p(alpa) = x_base + alpha x_dir,
  // with a bounding line.

  // Get parameterized line.
  Real x_base[2];
  for (int i = 0; i < 2; ++i)
    x_base[i] = 0.5*b/a[i];
  Real x_dir[] = {-a[1], a[0]};

  // Get the 4 alpha values.
  struct Alpha {
    Real alpha;
    Int side;
    bool operator<= (const Alpha& a) const { return alpha <= a.alpha; }
  };
  Alpha alphas[4];
  auto set_alpha = [&] (const Real* xbd, const Int& idx, const Int& side) {
    alphas[side].alpha = (xbd[idx] - x_base[idx])/x_dir[idx];
    alphas[side].side = side;
  };
  set_alpha(xlo, 1, 0); // bottom
  set_alpha(xhi, 0, 1); // right
  set_alpha(xhi, 1, 2); // top
  set_alpha(xlo, 0, 3); // left

  // Sort alphas. The middle two bound the feasible interval.
  sort4(alphas);

  // Eval the middle two and record the better of the them.
  auto eval_xi = [&] (const Real& alpha, const Int& idx) {
    return x_base[idx] + alpha*x_dir[idx];
  };
  auto eval_obj = [&] (const Real& alpha) {
    Real obj = 0;
    for (Int i = 0; i < 2; ++i) {
      x[i] = eval_xi(alpha, i);
      obj += w[i]*cedr::util::square(y[i] - x[i]);
    }
    return obj;
  };
  const Int ai = eval_obj(alphas[1].alpha) <= eval_obj(alphas[2].alpha) ? 1 : 2;

  Int info = 1, clipidx = 0;
  const auto& aai = alphas[ai];
  switch (aai.side) {
  case 0: x[0] = eval_xi(aai.alpha, 0); x[1] = xlo[1]; clipidx = 0; break;
  case 1: x[0] = xhi[0]; x[1] = eval_xi(aai.alpha, 1); clipidx = 1; break;
  case 2: x[0] = eval_xi(aai.alpha, 0); x[1] = xhi[1]; clipidx = 0; break;
  case 3: x[0] = xlo[0]; x[1] = eval_xi(aai.alpha, 1); clipidx = 1; break;
  default: cedr_assert(0); info = -2;
  }
  x[clipidx] = cedr::impl::min(xhi[clipidx], cedr::impl::max(xlo[clipidx], x[clipidx]));
  return info;
}

KOKKOS_INLINE_FUNCTION
Real calc_r_tol (const Real b, const Real* a, const Real* y, const Int n) {
  Real ab = std::abs(b);
  for (Int i = 0; i < n; ++i) ab = std::max(ab, std::abs(a[i]*y[i]));
  return 1e1*std::numeric_limits<Real>::epsilon()*std::abs(ab);
}

KOKKOS_INLINE_FUNCTION
void calc_r (const Int n, const Real* w, const Real* a, const Real b,
             const Real* xlo, const Real* xhi,  const Real* y, const Real& lambda,
             Real* x, Real& r, Real& r_lambda) {
  r = 0;
  r_lambda = 0;
  for (Int i = 0; i < n; ++i) {
    const Real q = a[i]/w[i];
    const Real x_trial = y[i] + lambda*q;
    Real xtmp;
    if (x_trial < (xtmp = xlo[i]))
      x[i] = xtmp;
    else if (x_trial > (xtmp = xhi[i]))
      x[i] = xtmp;
    else {
      x[i] = x_trial;
      r_lambda += a[i]*q;
    }
    r += a[i]*x[i];
  }
  r -= b;
}

KOKKOS_INLINE_FUNCTION
Int solve_1eq_bc_qp (const Int n, const Real* w, const Real* a, const Real b,
                     const Real* xlo, const Real* xhi, const Real* y, Real* x,
                     const Int max_its) {
  const Real r_tol = calc_r_tol(b, a, y, n);
  Int info;

  { // Check for a quick exit.
    bool all_in = true;
    Real r = 0;
    info = 0;
    for (Int i = 0; i < n; ++i) {
      if (x[i] != y[i]) {
        x[i] = y[i];
        info = 1;
      }
      if (is_outside(x[i], xlo, xhi, i)) {
        all_in = false;
        break;
      }
      r += a[i]*x[i];
    }
    if (all_in) {
      r -= b;
      if (std::abs(r) <= r_tol)
        return info;
    }
  }

  if (n == 2)
    return solve_1eq_bc_qp_2d(n, w, a, b, xlo, xhi, y, x);

  { // Eval r at end points to check for feasibility, and also possibly a quick
    // exit on a common case.
    Real r = -b;
    for (Int i = 0; i < n; ++i) {
      x[i] = xlo[i];
      r += a[i]*x[i];
    }
    if (std::abs(r) <= r_tol) return 1;
    if (r > 0) return -1;
    r = -b;
    for (Int i = 0; i < n; ++i) {
      x[i] = xhi[i];
      r += a[i]*x[i];
    }
    if (std::abs(r) <= r_tol) return 1;
    if (r < 0) return -1;
  }

  { // Check for a quick exit: the bounds are so tight that the midpoint of the
    // box satisfies r_tol.
    Real r = -b;
    for (Int i = 0; i < n; ++i) {
      x[i] = 0.5*(xlo[i] + xhi[i]);
      r += a[i]*x[i];
    }
    if (std::abs(r) <= r_tol) return 1;
  }

  const Real wall_dist = 1e-3;

  // Get lambda endpoints.
  Real lamlo = 0, lamhi = 0;
  for (Int i = 0; i < n; ++i) {
    const Real rq = w[i]/a[i];
    const Real lamlo_i = rq*(xlo[i] - y[i]);
    const Real lamhi_i = rq*(xhi[i] - y[i]);
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

  // Bisection-safeguarded Newton iteration for r(lambda) = 0.
  bool prev_step_bisect = false;
  Int nbisect = 0;
  info = -2;
  for (Int iteration = 0; iteration < max_its; ++iteration) {
    // Compute x, r, r_lambda.
    Real r, r_lambda;
    calc_r(n, w, a, b, xlo, xhi, y, lambda, x, r, r_lambda);
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

} // namespace local
} // namespace cedr

#endif
