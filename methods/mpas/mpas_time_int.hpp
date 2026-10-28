#ifndef INCLUDE_MPAS_TIME_INT_HPP
#define INCLUDE_MPAS_TIME_INT_HPP

#include "mpas_types.hpp"
#include "mpas_util.hpp"

#include <cstring>

namespace timeint {
struct OdeFn {
  virtual ~OdeFn () {}
  virtual bool eval(Real t, const Real* y, Real* f) const = 0;
  virtual void record (Real t, const Real* y) const {}
};

namespace simple {
void rk3(const Real ts, const Real dt, const int n, const Real* ys, Real* yp,
         const OdeFn& f, Real* ye);
// wrk has size >= 4*n.
void rk4(const Real ts, const Real dt, const int n, const Real* ys, Real* wrk,
         const OdeFn& f, Real* ye);
} // namespace simple

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

struct Info {
  Real good_initial_step;
};

struct ReturnState {
  enum Enum { success, function_eval_failed, step_too_small };
};

/*! \brief Implements the same RK3(2) pair as Matlab's ode23.
 *
 * A OdeFn f has
 * - method
 *       bool eval(Real t, const Real* y, Real* f) const
 *   to evaluate f(t), the ODE at time t. Return false on failure.
 * - method
 *       record(Real t, const Real* y)
 *   to optionally record y(t).
 *
 * \param opts [in] Options struct.
 * \param fun [in] ODE OdeFn.
 * \param y_caller [in/out] On input, y(t_s); on output, y(t_f).
 * \param n [in] length(y).
 * \param ts [in] t_s.
 * \param tf [in] t_f.
 * \param w [in/out] Workspace, of at least size 5n.
 *
 * Cite:
 *   P. Bogacki, L.F. Shampine, "A 3(2) Pair of Runge-Kutta Formulas",
 *   Appl. Math Lett. 2(4), 321-325, 1989.
 * and
 *   The MATLAB ODE Suite, L. F. Shampine and M. W. Reichelt, SIAM Journal on
 *   Scientific Computing, 18-1, 1997.
 */
ReturnState::Enum
ark23(const Options& opts, const OdeFn& fun,
      Real* const y_caller, const Size n,
      const Real ts, const Real tf,
      Real* const w, Info* info=0);

/*! \brief Implements the same RK5(4) pair as Matlab's ode45.
 *
 * \param w [in/out] Workspace, of at least size 9n.
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
ReturnState::Enum
ark45(const Options& opts, const OdeFn& fun,
      Real* const y_caller, const Size n,
      const Real ts, const Real tf,
      Real* const w, Info* info=0);

namespace test {
Int test_ark(const bool verbose);
} // namespace test
} // namespace timeint

#endif
