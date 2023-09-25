#ifndef INCLUDE_SLMM_TIME_INT_HPP
#define INCLUDE_SLMM_TIME_INT_HPP

#include "slmm_defs.hpp"
#include "slmm_util.hpp"

#include <cassert>
#include <cmath>
#include <cstring>
#include <limits>
#include <algorithm>
#include <iostream>
#include <iomanip>

namespace slmm {
namespace timeint {
class Options {
  Real initial_step_, rel_tol_, abs_tol_, max_step_size_;

public:
  Options ()
    : initial_step_(1e-3), rel_tol_(1e-3), abs_tol_(1e-6), max_step_size_(1e300)
  {}

  void set_initial_step (const Real is) { initial_step_ = is; }
  void set_rel_tol (const Real rt) { rel_tol_ = rt; }
  void set_abs_tol (const Real at) { abs_tol_ = at; }
  void set_max_step_size (const Real mss) { max_step_size_ = mss; }

  Real initial_step () const { return initial_step_; }
  Real rel_tol () const { return rel_tol_; }
  Real abs_tol () const { return abs_tol_; }
  Real max_step_size () const { return max_step_size_; }
};

struct Workspace {
  std::vector<Real> r;
};

struct Info {
  Real good_initial_step;
};

struct ReturnState {
  enum Enum { success, function_eval_failed, step_too_small };
};

template <typename T>
inline void copy (const Size n, const T* const s, T* const d)
{ for (Size i = 0; i < n; ++i) d[i] = s[i]; }

inline void aixiy (const Size n,
                   const Real a0, const Real* const x0,
                   const Real a1, const Real* const x1,
                   Real* const y) {
  for (Size i = 0; i < n; ++i)
    y[i] = a0*x0[i] + a1*x1[i];
}
inline void aixiy (const Size n,
                   const Real a0, const Real* const x0,
                   const Real a1, const Real* const x1,
                   const Real a2, const Real* const x2,
                   Real* const y) {
  for (Size i = 0; i < n; ++i)
    y[i] = a0*x0[i] + a1*x1[i] + a2*x2[i];
}
inline void aixiy (const Size n,
                   const Real a0, const Real* const x0,
                   const Real a1, const Real* const x1,
                   const Real a2, const Real* const x2,
                   const Real a3, const Real* const x3,
                   Real* const y) {
  for (Size i = 0; i < n; ++i)
    y[i] = a0*x0[i] + a1*x1[i] + a2*x2[i] + a3*x3[i];
}
inline void aixiy (const Size n,
                   const Real a0, const Real* const x0,
                   const Real a1, const Real* const x1,
                   const Real a2, const Real* const x2,
                   const Real a3, const Real* const x3,
                   const Real a4, const Real* const x4,
                   Real* const y) {
  for (Size i = 0; i < n; ++i)
    y[i] = a0*x0[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
}
inline void aixiy (const Size n,
                   const Real a0, const Real* const x0,
                   const Real a1, const Real* const x1,
                   const Real a2, const Real* const x2,
                   const Real a3, const Real* const x3,
                   const Real a4, const Real* const x4,
                   const Real a5, const Real* const x5,
                   Real* const y) {
  for (Size i = 0; i < n; ++i)
    y[i] = a0*x0[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i] + a5*x5[i];
}

/*! \brief Implements the same RK3(2) pair as Matlab's ode23.
 *
 * A Functor f has
 * - method
 *       bool eval(Real t, const Real* y, Real* f) const
 *   to evaluate f(t), the ODE at time t. Return false on failure.
 * - method
 *       record(Real t, const Real* y) const
 *   to optionally record y(t).
 *
 * \param opts [in] Options struct.
 * \param fun [in] ODE Functor.
 * \param y_caller [in/out] On input, y(t_s); on output, y(t_f).
 * \param n [in] length(y).
 * \param ts [in] t_s.
 * \param tf [in] t_f.
 * \param w [in/out] Workspace. Reuse between calls to minimize allocations.
 *
 * Cite:
 *   P. Bogacki, L.F. Shampine, "A 3(2) Pair of Runge-Kutta Formulas",
 *   Appl. Math Lett. 2(4), 321-325, 1989.
 * and
 *   The MATLAB ODE Suite, L. F. Shampine and M. W. Reichelt, SIAM Journal on
 *   Scientific Computing, 18-1, 1997.
 */
template<typename Functor>
ReturnState::Enum
ark23 (const Options& opts, const Functor& fun,
       Real* const y_caller, const Size n,
       const Real ts, const Real tf,
       Workspace& w, Info* info=0) {
  static const Real pow = 1.0/3.0;

  const Real threshold = opts.abs_tol() / opts.rel_tol();
  const int tdir = tf >= ts ? 1 : -1;

  w.r.resize(5*n);
  Real* f0 = w.r.data();
  Real* const f1 = f0 + n;
  Real* const f2 = f1 + n;
  Real* f3 = f2 + n;
  Real* y0 = y_caller;
  Real* y1 = f3 + n;

  Real t = ts;
  const Real sgn = sign(tf - ts);
  Real absh = std::abs(opts.initial_step());
  if (info) info->good_initial_step = absh;
  bool fgood = fun.eval(t, y0, f0);
  fun.record(t, y0);
  if ( ! fgood) return ReturnState::function_eval_failed;

  while (sgn*t < sgn*tf) {
    const double hmin = 16*std::numeric_limits<Real>::epsilon()*t;
    bool no_failed = true;
    Real err, tnew;
    for (;;) { // Integrate one step; loop until success.
      // Get tnew and sanitized h.
      absh = std::min(absh, opts.max_step_size());
      Real h = tdir*absh;
      if (sgn*(t + h) > sgn*tf) {
        h = tf - t;
        absh = std::abs(h);
      }
      tnew = t + h;
      h = tnew - t;

      // Integration rule.
      do {
        aixiy(n, 1, y0, 0.5*h, f0, y1);
        fgood = fun.eval(t + 0.5*h, y1, f1);
        if ( ! fgood) break;
        aixiy(n, 1, y0, 0.75*h, f1, y1);
        fgood = fun.eval(t + 0.75*h, y1, f2);
        if ( ! fgood) break;
        aixiy(n, 1, y0, 2.0*h/9.0, f0, h/3.0, f1, 4.0*h/9.0, f2, y1);
        fgood = fun.eval(tnew, y1, f3);
      } while (0);

      // Determine error.
      err = 0;
      if ( ! fgood) {
        err = opts.rel_tol() + 1;
        no_failed = false;
      } else {
        // Coefficients from subtracting the order-2 prediction from the order-3
        // prediction:
        static const Real E[] = {-5.0/72.0, 1.0/12.0, 1.0/9.0, -1.0/8.0};
        // Element-wise error control:
        //   err = absh * norm( (f*E) ./ max( max(abs(y), abs(yt)),
        //                                    threshold),
        //                      inf );
        for (Size i = 0; i < n; ++i) {
          const Real fE =
            std::abs(E[0]*f0[i] + E[1]*f1[i] + E[2]*f2[i] + E[3]*f3[i]);
          const Real den =
            std::max<Real>(std::max<Real>(std::abs(y0[i]), std::abs(y1[i])),
                           threshold);
          err = std::max<Real>(err, fE / den);
        }
        err *= absh;
      }

      // Determine if the step succeeded. If it did not, compute a smaller step
      // size and try again.
      if (err > opts.rel_tol()) {
        if (absh <= hmin) {
          fun.record(t, y1);
          return ReturnState::step_too_small;
        }
        if (no_failed) {
          no_failed = false;
          absh = std::max(
            hmin, absh*std::max<Real>(
              0.5, 0.8*std::pow(opts.rel_tol()/err, pow)));
        } else {
          absh = std::max<Real>(hmin, 0.5*absh);
        }
      } else {
        // Successful step. Break from the integration loop.
        break;
      }
    } // One integration step.
    if (info) info->good_initial_step = absh;

    if (no_failed) {
      // Integration step succeeded on first try. Increase the step size.
      const Real fac = 0.8*std::pow(opts.rel_tol()/err, pow);
      // Don't increase the step size by more than 5x.
      absh = std::min<Real>(5, fac)*absh;
    }

    t = tnew;
    // Swap pointers.
    std::swap(y0, y1);
    std::swap(f0, f3);

    fun.record(t, y0);
  }

  // On output, y_caller contains y(tf). If the pointers don't agree (because of
  // swapping above), copy.
  if (y_caller != y0)
    memcpy(y_caller, y0, n*sizeof(*y_caller));

  return ReturnState::success;
}

/*! \brief Implements the same RK5(4) pair as Matlab's ode45.
 *
 * Cite:
 *   Dormand, J. R.; Prince, P. J. (1980), "A family of embedded Runge-Kutta
 *   formulae", Journal of Computational and Applied Mathematics 6 (1): 19–26.
 * and
 *   The MATLAB ODE Suite, L. F. Shampine and M. W. Reichelt, SIAM Journal on
 *   Scientific Computing, 18-1, 1997.
 *
 * The Butcher tableau is
 *
 * 0    |
 * 1/5  | 1/5
 * 3/10 | 3/40         9/40
 * 4/5  | 44/45       -56/15        32/9
 * 8/9  | 19372/656   -25360/2187   64448/6561   -212/729
 * 1    | 9017/3168   -355/33       46732/5247    49/176    -5103/18656
 * 1    | 35/384       0            500/1113      125/192   -2187/6784      11/84   
 * --------------------------------------------------------------------------------
 *      | 35/384       0            500/1113      125/192   -2187/6784      11/84     0
 *      | 5179/57600   0            7571/16695    393/640   -92097/339200   187/2100  1/40
 *
 * and the corresponding E array, obtained from subtracting the first row of b
 * from the second, is
 *
 *       -71/57600     0            71/16695     -71/1920    17253/339200  -88/2100   1/40
 */
template<typename Functor>
ReturnState::Enum
ark45 (const Options& opts, const Functor& fun,
       Real* const y_caller, const Size n,
       const Real ts, const Real tf,
       Workspace& w, Info* info=0) {
  static const Real
    c2 = 0.2, c3 = 0.3, c4 = 0.8, c5 = 8.0/9.0;
  static const Real
    a21 = c2,
    a31 = 3.0/40.0, a32 = 9.0/40.0,
    a41 = 44.0/45.0, a42 = -56.0/15.0, a43 = 32.0/9.0,
    a51 = 19372.0/6561.0, a52 = -25360.0/2187.0, a53 = 64448.0/6561.0,
      a54 = -212.0/729.0, 
    a61 = 9017.0/3168.0, a62 = -355.0/33.0, a63 = 46732.0/5247.0,
      a64 = 49.0/176.0, a65 = -5103.0/18656.0, 
    a71 = 35.0/384.0, a73 = 500.0/1113.0, a74 = 125.0/192.0,
      a75 = -2187.0/6784.0, a76 = 11.0/84.0;
  static const Real pow = 1.0/5.0;
  // Coefficients from subtracting the order-4 prediction from the order-5
  // prediction:
  static const Real E[] = {-71.0/57600.0, 0.0, 71.0/16695.0, -71.0/1920.0,
                           17253.0/339200.0, -88.0/2100.0, 1.0/40.0};

  const Real threshold = opts.abs_tol() / opts.rel_tol();
  const int tdir = tf >= ts ? 1 : -1;

  w.r.resize(9*n);
  Real*       f0 = w.r.data();
  Real* const f1 = f0 + n;
  Real* const f2 = f1 + n;
  Real* const f3 = f2 + n;
  Real* const f4 = f3 + n;
  Real* const f5 = f4 + n;
  Real*       f6 = f5 + n;
  Real*       y0 = y_caller;
  Real*       y1 = f6 + n;

  Real t = ts;
  const Real sgn = sign(tf - ts);
  Real absh = std::abs(opts.initial_step());
  if (info) info->good_initial_step = absh;
  bool fgood = fun.eval(t, y0, f0);
  fun.record(t, y0);
  if ( ! fgood) return ReturnState::function_eval_failed;

  while (sgn*t < sgn*tf) {
    const double hmin = 16*std::numeric_limits<Real>::epsilon()*t;
    bool no_failed = true;
    Real err, tnew;
    for (;;) { // Integrate one step; loop until success.
      // Get tnew and sanitized h.
      absh = std::min(absh, opts.max_step_size());
      Real h = tdir*absh;
      if (sgn*(t + h) > sgn*tf) {
        h = tf - t;
        absh = std::abs(h);
      }
      tnew = t + h;
      h = tnew - t;

      // Integration rule.
      do {
        aixiy(n, 1, y0, a21*h, f0, y1);
        fgood = fun.eval(t + c2*h, y1, f1);
        if ( ! fgood) break;
        aixiy(n, 1, y0, a31*h, f0, a32*h, f1, y1);
        fgood = fun.eval(t + c3*h, y1, f2);
        if ( ! fgood) break;
        aixiy(n, 1, y0, a41*h, f0, a42*h, f1, a43*h, f2, y1);
        fgood = fun.eval(t + c4*h, y1, f3);
        if ( ! fgood) break;
        aixiy(n, 1, y0, a51*h, f0, a52*h, f1, a53*h, f2, a54*h, f3, y1);
        fgood = fun.eval(t + c5*h, y1, f4);
        if ( ! fgood) break;
        aixiy(n, 1, y0, a61*h, f0, a62*h, f1, a63*h, f2, a64*h, f3, a65*h, f4, y1);
        fgood = fun.eval(tnew,     y1, f5);
        if ( ! fgood) break;
        aixiy(n, 1, y0, a71*h, f0,            a73*h, f2, a74*h, f3, a75*h, f4, a76*h, f5, y1);
        fgood = fun.eval(tnew,     y1, f6);
      } while (0);

      // Determine error.
      err = 0;
      if ( ! fgood) {
        err = opts.rel_tol() + 1;
        no_failed = false;
      } else {
        for (Size i = 0; i < n; ++i) {
          const Real fE =
            std::abs(E[0]*f0[i] + E[1]*f1[i] + E[2]*f2[i] + E[3]*f3[i] +
                     E[4]*f4[i] + E[5]*f5[i] + E[6]*f6[i]);
          const Real den =
            std::max<Real>(std::max<Real>(std::abs(y0[i]), std::abs(y1[i])),
                           threshold);
          err = std::max<Real>(err, fE / den);
        }
        err *= absh;
      }

      // Determine if the step succeeded. If it did not, compute a smaller step
      // size and try again.
      if (err > opts.rel_tol()) {
        if (absh <= hmin) {
          fun.record(t, y1);
          return ReturnState::step_too_small;
        }
        if (no_failed) {
          no_failed = false;
          absh = std::max(
            hmin, absh*std::max<Real>(
              0.5, 0.8*std::pow(opts.rel_tol()/err, pow)));
        } else {
          absh = std::max<Real>(hmin, 0.5*absh);
        }
      } else {
        // Successful step. Break from the integration loop.
        break;
      }
    } // One integration step.
    if (info) info->good_initial_step = absh;

    if (no_failed) {
      // Integration step succeeded on first try. Increase the step size.
      const Real fac = 0.8*std::pow(opts.rel_tol()/err, pow);
      // Don't increase the step size by more than 5x.
      absh = std::min<Real>(5, fac)*absh;
    }

    t = tnew;
    // Swap pointers.
    std::swap(y0, y1);
    std::swap(f0, f6);

    fun.record(t, y0);
  }

  // On output, y_caller contains y(tf). If the pointers don't agree (because of
  // swapping above), copy.
  if (y_caller != y0)
    memcpy(y_caller, y0, n*sizeof(*y_caller));

  return ReturnState::success;
}

namespace test {
Int test_ark(const bool verbose);
} // namespace test
} // namespace timeint
} // namespace slmm

#endif
