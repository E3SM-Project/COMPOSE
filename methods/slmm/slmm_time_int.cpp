#include "slmm_time_int.hpp"
#include "slmm_util.hpp"

namespace slmm {
namespace timeint {
namespace test {
class TestFunctor {
  mutable Size nsteps_;

protected:
  Real tspan_[2], ys_[2];

public:
  TestFunctor (const Real tspan[2], const Real ys[2])
    : nsteps_(0)
  {
    copy(2, tspan, tspan_);
    copy(2, ys, ys_);
  }
  Size nsteps () const { return nsteps_; }
  void reset () { nsteps_ = 0; }
  const Real* tspan () const { return tspan_; }
  const Real* ys () const { return ys_; }
  void record (const Real t, const Real* const y) const { ++nsteps_; }
  virtual bool eval (const Real t, const Real* const y, Real* const f) const = 0;
  virtual void eval_solution(const Real t, Real* const f) const = 0;
};

// ODE
//     y'(t) = lambda y(t)
// with solution
//     y(tf) = y(ts) e^(lambda (tf - ts)).
class LambdaFunctor : public TestFunctor {
  Real lambda_[2];
public:
  LambdaFunctor (const Real lambda[2], const Real tspan[2], const Real ys[2])
    : TestFunctor(tspan, ys)
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

class TimeDepFunctor : public TestFunctor {
  const Real a_;
public:
  TimeDepFunctor (const Real a, const Real tspan[2], const Real ys[2])
    : TestFunctor(tspan, ys), a_(a)
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

bool test_ark_y2 (TestFunctor& fun, const ARKMethod method,
                  const bool verbose = true) {
  std::ostream& os = std::cout;
  auto ios_state = save_ios(os);

  Options opts;
  opts.set_initial_step(1e-3);
  opts.set_abs_tol(1e-20);
  Workspace w;
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
      ark23(opts, fun, y, 2, fun.tspan()[0], fun.tspan()[1], w);
    else
      ark45(opts, fun, y, 2, fun.tspan()[0], fun.tspan()[1], w);

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
      LambdaFunctor fun(lambda, tspan, ys);
      pass = pass && test_ark_y2(fun, method_ark23, verbose);
      pass = pass && test_ark_y2(fun, method_ark45, verbose);
    }
    {
      if (verbose) std::cout << "  Standard test function backwards in time.\n";
      const Real tspanb[] = {8, -3};
      LambdaFunctor fun(lambda, tspanb, ys);
      pass = pass && test_ark_y2(fun, method_ark23, verbose);
      pass = pass && test_ark_y2(fun, method_ark45, verbose);
    }
  }
  {
    if (verbose) std::cout << "  Exact time-dependent function.\n";
    TimeDepFunctor fun(0.1, tspan, ys);
    pass = pass && test_ark_y2(fun, method_ark23, verbose);
    pass = pass && test_ark_y2(fun, method_ark45, verbose);
  }
  return pass ? 0 : 1;
}

} // namespace test
} // namespace timeint
} // namespace slmm
