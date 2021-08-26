#include "islet_np4.hpp"
#include "islet_isl.hpp"

static const Real oosqrt5 = 0.44721359549995793928;

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

Np4InterpMethod::Np4InterpMethod (Real c0, Real c1, Real c2) {
  reset_c(c0, c1, c2);
}

void Np4InterpMethod::reset_c (Real c0, Real c1, Real c2) {
  c[0] = c0; c[1] = c1; c[2] = c2;
}

Real Np4InterpMethod::eval_a (const Real& x) const {
  return eval_lagrange_poly(3, islet::get_x_gll(3), c,
                            2*(1 - std::abs(x))/(1 - oosqrt5) - 1);
}

void Np4InterpMethod::eval (const Real& x, Real* const y) {
  const auto* x_gll = get_xnodes();
  if (x < -oosqrt5 || x > oosqrt5) {
    y[0] = y[3] = 0;
    const Int os = x < -oosqrt5 ? 0 : 1;
    islet::eval_lagrange_poly(x_gll + os, 3, x, y + os);
    Real y4[4];
    islet::eval_lagrange_poly(x_gll, 4, x, y4);
    const Real a = eval_a(x);
    for (int i = 0; i < 4; ++i)
      y[i] = a*y[i] + (1 - a)*y4[i];
  } else
    islet::eval_lagrange_poly(x_gll, 4, x, y);
}
