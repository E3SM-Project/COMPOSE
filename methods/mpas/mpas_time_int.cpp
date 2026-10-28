#include "mpas_time_int.hpp"

#include <iostream>
#include <iomanip>

template<typename Stream> class IosSaver {
  Stream& s_;
  std::ios state_;
public:
  IosSaver (Stream& s) : s_(s), state_(nullptr) { state_.copyfmt(s); }
  IosSaver (const IosSaver& ios) : s_(ios.s_), state_(nullptr)
  { state_.copyfmt(ios.state_); }
  IosSaver operator= (const IosSaver&) = delete;
  ~IosSaver () { s_.copyfmt(state_); }
};
template<typename Stream> inline IosSaver<Stream> save_ios (Stream& s)
{ return IosSaver<Stream>(s); }

namespace timeint {
namespace simple {
// Eq 2 in Wang, Skamarock, Feingold, MWR, 2009 and eq 3 in Skamarock and
// Gassmann, MWR, 2011.
void rk3 (const Real ts, const Real dt, const int n, const Real* ys, Real* yp,
          const OdeFn& f, Real* ye) {
  static const Real dtcoef[] = {0, 1.0/3.0, 0.5, 1};
  for (int it = 0; it < 3; ++it) {
    f.eval(ts + dtcoef[it]*dt, it == 0 ? ys : ye, yp);
    for (int i = 0; i < n; ++i) ye[i] = ys[i] + dtcoef[it+1]*dt*yp[i];
  }
}

void rk4 (const Real ts, const Real dt, const int n, const Real* ys, Real* wrk,
          const OdeFn& f, Real* ye) {
  static const Real dtcoef[] = {0, 0.5, 0.5, 1, 1};
  static const Real b[] = {1.0/6.0, 2.0/6.0, 2.0/6.0, 1.0/6.0};
  Real* const yp = wrk;
  wrk += n;
  for (int it = 0; it < 4; ++it) {
    f.eval(ts + dtcoef[it]*dt, it == 0 ? ys : ye, yp);
    if (it < 3) {
      for (int i = 0; i < n; ++i)
        ye[i] = ys[i] + dtcoef[it+1]*dt*yp[i];
      for (int i = 0; i < n; ++i)
        wrk[n*it+i] = yp[i];
    } else {
      for (int i = 0; i < n; ++i)
        ye[i] = ys[i] + dt*(b[0]*wrk[i] + b[1]*wrk[n+i] + b[2]*wrk[2*n+i] + b[3]*yp[i]);
    }
  }
}
} // namespace simple

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

ReturnState::Enum
ark23 (const Options& opts, const OdeFn& fun,
       Real* const y_caller, const Size n,
       const Real ts, const Real tf,
       Real* const w, Info* info) {
  static const Real pow = 1.0/3.0;

  const Real threshold = opts.abs_tol() / opts.rel_tol();
  const int tdir = tf >= ts ? 1 : -1;

  Real* f0 = w;
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
      const Real fac = err > 0 ? 0.8*std::pow(opts.rel_tol()/err, pow) : 5;
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

ReturnState::Enum
ark45 (const Options& opts, const OdeFn& fun,
       Real* const y_caller, const Size n,
       const Real ts, const Real tf,
       Real* const w, Info* info) {
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

  Real*       f0 = w;
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
      const Real fac = err > 0 ? 0.8*std::pow(opts.rel_tol()/err, pow) : 5;
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
class TestOdeFn : public OdeFn {
  mutable Size nsteps_;

protected:
  Real tspan_[2], ys_[2];

public:
  TestOdeFn (const Real tspan[2], const Real ys[2])
    : nsteps_(0)
  {
    copy(2, tspan, tspan_);
    copy(2, ys, ys_);
  }
  Size nsteps () const { return nsteps_; }
  void reset () { nsteps_ = 0; }
  const Real* tspan () const { return tspan_; }
  const Real* ys () const { return ys_; }
  void record (const Real t, const Real* const y) const override { ++nsteps_; }
  virtual void eval_solution(const Real t, Real* const f) const = 0;
};

// ODE
//     y'(t) = lambda y(t)
// with solution
//     y(tf) = y(ts) e^(lambda (tf - ts)).
class LambdaOdeFn : public TestOdeFn {
  Real lambda_[2];
public:
  LambdaOdeFn (const Real lambda[2], const Real tspan[2], const Real ys[2])
    : TestOdeFn(tspan, ys)
  {
    copy(2, lambda, lambda_);
  }
  virtual bool eval (const Real t, const Real* const y, Real* const f) const {
    f[0] = lambda_[0]*y[0] - lambda_[1]*y[1];
    f[1] = lambda_[0]*y[1] + lambda_[1]*y[0];
    return true;
  }
  virtual void eval_solution (const Real t, Real* const y) const {
    const Real dt = t - tspan_[0];
    const Real
      c = std::cos(dt*lambda_[1]),
      s = std::sin(dt*lambda_[1]),
      mag = std::exp(dt*lambda_[0]);
    y[0] = mag*(ys_[0]*c - ys_[1]*s);
    y[1] = mag*(ys_[0]*s + ys_[1]*c);
  }
};

class TimeDepOdeFn : public TestOdeFn {
  const Real a_;
public:
  TimeDepOdeFn (const Real a, const Real tspan[2], const Real ys[2])
    : TestOdeFn(tspan, ys), a_(a)
  {}
  virtual bool eval (const Real t, const Real* const y, Real* const f) const {
    f[0] = a_*t;
    f[1] = -0.5*a_*t;
    return true;
  }
  virtual void eval_solution (const Real t, Real* const y) const {
    const Real dst = square(t) - square(tspan_[0]);
    y[0] = ys_[0] +  0.5*a_*dst;
    y[1] = ys_[1] - 0.25*a_*dst;
  }
};

enum ARKMethod { method_ark23, method_ark45 };

bool test_ark_y2 (TestOdeFn& fun, const ARKMethod method,
                  const bool verbose = true) {
  std::ostream& os = std::cout;
  auto ios_state = save_ios(os);

  Options opts;
  opts.set_initial_step(1e-3);
  opts.set_abs_tol(1e-20);
  Real wrk[18];
  Real ya[2];
  fun.eval_solution(fun.tspan()[1], ya);

  Real rtol = 1e-1;
  Real rds[6];
  const Size ntrial = static_cast<Size>(sizeof(rds)/sizeof(*rds));
  const Real rtol_increase = 100;
  const Real den = std::sqrt(square(ya[0]) + square(ya[1]));

  for (Size trial = 0; trial < ntrial; ++trial) {
    rtol *= 1/rtol_increase;
    opts.set_rel_tol(rtol);
    fun.reset();
    Real y[2];
    copy(2, fun.ys(), y);
    if (method == method_ark23)
      ark23(opts, fun, y, 2, fun.tspan()[0], fun.tspan()[1], wrk);
    else
      ark45(opts, fun, y, 2, fun.tspan()[0], fun.tspan()[1], wrk);

    rds[trial] = std::sqrt(square(y[0] - ya[0]) + square(y[1] - ya[1])) / den;
    if (verbose) {
      os.precision(2);
      os << "  trial " << std::setw(2) << trial
         << " nsteps " << std::setw(6) << fun.nsteps()
         << " rtol " << std::scientific << rtol
         << " reldif " << rds[trial] << "\n";
    }
  }

  const Real improvement = rds[ntrial-2] / rds[ntrial-1];
  const bool pass = 
    rds[ntrial-1] <= 1e2*rtol &&
    (improvement >= 0.9*rtol_increase ||
     rds[ntrial-1] <= 1e2*std::numeric_limits<Real>::epsilon());
  return pass;
}

Int test_ark (const bool verbose) {
  if (verbose)
    std::cout << "> Adaptive Runge-Kutta 2-3 unit test\n";
  static const Real tspan[] = {0.5, 71.2}, ys[] = {3.6, -0.7};
  bool pass = true;
  {
    static const Real lambda[] = {-0.02, 0.25};
    {
      if (verbose) std::cout << "  Standard test function.\n";
      LambdaOdeFn fun(lambda, tspan, ys);
      pass = pass && test_ark_y2(fun, method_ark23, verbose);
      pass = pass && test_ark_y2(fun, method_ark45, verbose);
    }
    {
      if (verbose) std::cout << "  Standard test function backwards in time.\n";
      const Real tspanb[] = {8, -3};
      LambdaOdeFn fun(lambda, tspanb, ys);
      pass = pass && test_ark_y2(fun, method_ark23, verbose);
      pass = pass && test_ark_y2(fun, method_ark45, verbose);
    }
  }
  {
    if (verbose) std::cout << "  Exact time-dependent function.\n";
    TimeDepOdeFn fun(0.1, tspan, ys);
    pass = pass && test_ark_y2(fun, method_ark23, verbose);
    pass = pass && test_ark_y2(fun, method_ark45, verbose);
  }
  return pass ? 0 : 1;
}

} // namespace test
} // namespace timeint
