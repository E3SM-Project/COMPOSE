#include "slmm_islet.hpp"

namespace slmm {
namespace islet {

bool GllOffsetNodal::get_w (const Int& np, const Real*& w) const {
  switch (np) {
  case 2: w = w_np2; break;
  case 3: w = w_np3; break;
  case 4: w = w_np4; break;
  case 5: w = w_np5; break;
  case 6: w = w_np6; break;
  case 7: w = w_np7; break;
  case 8: w = w_np8; break;
  case 9: w = w_np9; break;
  case 10: w = w_np10; break;
  case 11: w = w_np11; break;
  case 12: w = w_np12; break;
  case 13: w = w_np13; break;
  case 16: w = w_np16; break;
  default: w = nullptr; return false;
  }
  return true;
}

Int GllOffsetNodal::max_degree (const Int& np) const {
  static const Int degrees[] = {-1,-1,1,2,3,3,5,5,6,7,7,8,9,10,-1,-1,12};
  return degrees[np];
}

static Real normalize_x (const Real* gll_x, const Real& x) {
  const Real x0 = gll_x[1];
  return (x - x0) / (1 - x0);
}

static void outer_eval (const Real* gll_x, const Real& x, Real v[4]) {
  const Real
    xbar = normalize_x(gll_x, gll_x[2]),
    ooxbar = 1 / xbar,
    ybar = 1 / (xbar - 1),
    xn = normalize_x(gll_x, x);
  v[0] = 0;
  v[1] = 1 + ybar*xn*((1 - ooxbar)*xn + ooxbar - xbar);
  v[2] = ybar*ooxbar*xn*(xn - 1);
  v[3] = ybar*xn*(xbar - xn);
}

#if 0
static bool np4_subgrid_eval (const Real* const x_gll, const Real& x,
                              Real y[4]) {
  static constexpr Real
    alpha = 0.5527864045000416708,
    v = 0.427*(1 + alpha),
    x2 = 0.4472135954999579277, // 1/sqrt(5)
    x3 = 1 - x2,
    det = x2*x3*(x2 - x3),
    y2 = alpha,
    y3 = v,
    c1 = (x3*y2 - x2*y3)/det,
    c2 = (-x3*x3*y2 + x2*x2*y3)/det;
  if (x < x_gll[1] || x > x_gll[2]) {
    Real y4[4];
    GLL::eval_lagrange_poly(4, x_gll, x, y4);
    if (x < x_gll[1]) {
      outer_eval(x_gll, -x, y);
      std::swap(y[0], y[3]);
      std::swap(y[1], y[2]);
    } else
      outer_eval(x_gll, x, y);
    const Real x0 = 1 - std::abs(x);
    const Real a = (c1*x0 + c2)*x0;
    for (int i = 0; i < 4; ++i)
      y[i] = a*y[i] + (1 - a)*y4[i];
  }
  else
    GLL::eval_lagrange_poly(4, x_gll, x, y);
  return true;
}
#else
static bool np4_subgrid_eval (const Real* const x_gll, const Real& x,
                              Real y[4]) {
  static const Real c1 = 0.306;
  if (x < x_gll[1] || x > x_gll[2]) {
    y[0] = y[3] = 0;
    const Int os = x < x_gll[1] ? 0 : 1;
    GLL::eval_lagrange_poly(3, x_gll + os, x, y + os);
    Real y4[4];
    GLL::eval_lagrange_poly(4, x_gll, x, y4);
    const Real x0 = 2*(1 - std::abs(x))/(1 - x_gll[2]) - 1;
    const Real a = (c1 + (0.5 - c1)*x0)*(x0 + 1);
    for (int i = 0; i < 4; ++i)
      y[i] = a*y[i] + (1 - a)*y4[i];
  } else
    GLL::eval_lagrange_poly(4, x_gll, x, y);
  return true;
}
#endif

bool GllOffsetNodal::eval (const Int& np, const Real& x, Real* const v) const {
  const Real* xnode;
  get_x(np, xnode);
  switch (np) {
  case  2: return evalon< 2,0>(xnode, {                      }, {                }, x, v); // 1
  case  3: return evalon< 3,0>(xnode, {                      }, {                }, x, v); // 2
  case  4: return np4_subgrid_eval(xnode, x, v);                                           // 2
  case  5: return evalon< 5,2>(xnode, { 3,  4                }, {0, 0            }, x, v); // 2
  case  6: return evalon< 6,2>(xnode, { 6,  5                }, {0, 0            }, x, v); // 4
  case  7: return evalon< 7,3>(xnode, { 5,  5,  6            }, {0, 0, 0         }, x, v); // 4
  case  8: return evalon< 8,4>(xnode, { 6,  6,  7,  6        }, {0, 0, 0, 1      }, x, v); // 5
  case  9: return evalon< 9,4>(xnode, { 7,  8,  7,  7        }, {0, 0, 0, 1      }, x, v); // 6
  case 10: return evalon<10,5>(xnode, { 7,  7,  7,  8,  8    }, {0, 0, 0, 0, 1   }, x, v); // 6
  case 11: return evalon<11,5>(xnode, { 8,  9,  8,  9,  8    }, {0, 0, 0, 0, 1   }, x, v); // 7
  case 12: return evalon<12,6>(xnode, { 9,  9, 10, 10,  9, 10}, {0, 0, 0, 0, 1, 1}, x, v); // 8
  case 13: return evalon<13,6>(xnode, {10, 10, 10, 10, 11, 10}, {0, 0, 0, 0, 0, 1}, x, v); // 9
  case 16: return evalon<16,8>(xnode, {12, 13, 13, 13, 13, 14, 13, 12},
                                      { 0,  0,  0,  0,  0,  0,  1,  2}, x, v); // 11
  default: return false;
  }
  return true;
}

bool GllOffsetNodal::eval_derivative (const Int& np, const Real& x, Real* const v) const {
  assert(0);
  return false;
}

bool GllNodal::get_w (const Int& np, const Real*& w) const {
  switch (np) {
  case 2: w = w_np2; break;
  case 3: w = w_np3; break;
  case 4: w = w_np4; break;
  case 5: w = w_np5; break;
  case 6: w = w_np6; break;
  case 7: w = w_np7; break;
  case 8: w = w_np8; break;
  case 9: w = w_np9; break;
  case 10: w = w_np10; break;
  case 11: w = w_np11; break;
  case 12: w = w_np12; break;
  case 13: w = w_np13; break;
  case 16: w = GllOffsetNodal::w_np16; break;
  default: w = nullptr; return false;
  }
  return true;
}

void eval (const Int& np, const Real* const xnodes,
           const Int* subnp, Int const* const* nodes,
           const Real& x, Real* const v) {
  if (x > 0) {
    eval(np, xnodes, subnp, nodes, -x, v);
    for (int i = 0; i < np/2; ++i)
      std::swap(v[i], v[np-i-1]);
    return;
  }
  Real xsub[Basis::np_max], vsub[Basis::np_max];
  const Int nreg = np-1, ios = 1;
  for (Int i = 0; i < nreg; ++i) {
    if (i < np-2 && x > xnodes[i+ios]) continue;
    if (subnp[i] == np) {
      Basis::eval_lagrange_poly(np, xnodes, x, v);
    } else {
      for (Int j = 0; j < subnp[i]; ++j)
        xsub[j] = xnodes[nodes[i][j]];
      std::fill(v, v + np, 0);
      Basis::eval_lagrange_poly(subnp[i], xsub, x, vsub);
      for (Int j = 0; j < subnp[i]; ++j) {
        const auto node = nodes[i][j];
        assert(node >= 0);
        assert(node < np);
        v[node] = vsub[j];
      }
    }
    break;
  }
}

bool GllNodal::eval (const Int& np, const Real& x, Real* const v) const {
  switch (np) {
  case 6: { // 4
    const Real* xnodes;
    get_x(np, xnodes);
    const Int subnp[] = {5,5,6};
    const Int n0[] = { 0, 1, 2, 3, 4,  };
    const Int n1[] = { 0, 1, 2, 3,    5};
    const Int n2[] = { 0, 1, 2, 3, 4, 5};
    const Int* nodes[] = {n0,n1,n2};
    islet::eval(6, xnodes, subnp, nodes, x, v);
  } break;
  case 9: { // 6
    const Real* xnodes;
    get_x(np, xnodes);
    const Int subnp[] = {7,8,8,7};
    const Int n0[] = { 0, 1, 2, 3, 4, 5,       8};
    const Int n1[] = { 0, 1, 2, 3, 4, 5,    7, 8};
    const Int n2[] = { 0, 1, 2, 3, 4, 5, 6,    8};
    const Int n3[] = {    1, 2, 3, 4, 5, 6, 7   };
    const Int* nodes[] = {n0,n1,n2,n3};
    islet::eval(9, xnodes, subnp, nodes, x, v);
  } break;
  default: return GllOffsetNodal::eval(np, x, v);
  }
  return true;
}

bool GllNodal::eval_derivative (const Int& np, const Real& x, Real* const v) const {
  assert(0);
  return false;
}

Int UniformOffsetNodal::max_degree (const Int& np) const { return 4; }

bool UniformOffsetNodal::get_x (const Int& np, const Real*& x) const {
  // Provide x data for some bases we don't actually support.
  if (np < 2 || np > Basis::np_max+1) {
    x = nullptr;
    return false;
  }
  static Real xnode[Basis::np_max+1][Basis::np_max+1] = {0};
  if (xnode[np][0] == 0) {
    for (Int i = 0; i < np; ++i)
      xnode[np][i] = 2*(Real(i)/(np-1)) - 1;
  }
  x = xnode[np];
  return true;
}

bool UniformOffsetNodal::get_w (const Int& np, const Real*& w) const {
  switch (np) {
  case 2: w = w_np2; break;
  case 3: w = w_np3; break;
  case 4: w = w_np4; break;
  case 5: w = w_np5; break;
  case 6: w = w_np6; break;
  case 7: w = w_np7; break;
  case 8: w = w_np8; break;
  case 9: w = w_np9; break;
  case 10: w = w_np10; break;
  case 11: w = w_np11; break;
  case 12: w = w_np12; break;
  case 13: w = w_np13; break;
  default: w = nullptr; return false;
  }
  return true;
}

bool UniformOffsetNodal
::eval (const Int& np, const Real& x, Real* const v) const {
  const Real* xnode;
  get_x(np, xnode);
  switch (np) {
  case  2: return evalon< 2,0>(xnode, {           }, {           }, x, v);
  case  3: return evalon< 3,0>(xnode, {           }, {           }, x, v);
  case  4: return evalon< 4,2>(xnode, {3,4        }, {0,0        }, x, v);
  case  5: return evalon< 5,2>(xnode, {3,4        }, {0,0        }, x, v);
  case  6: return evalon< 6,3>(xnode, {3,4,6      }, {0,0,0      }, x, v);
  case  7: return evalon< 7,3>(xnode, {3,4,4      }, {0,0,1      }, x, v);
  case  8: return evalon< 8,4>(xnode, {4,4,4,4    }, {0,0,1,2    }, x, v);
  case  9: return evalon< 9,4>(xnode, {4,4,4,4    }, {0,0,1,2    }, x, v);
  case 10: return evalon<10,5>(xnode, {4,4,4,4,4  }, {0,0,1,2,3  }, x, v);
  case 11: return evalon<11,5>(xnode, {4,4,4,4,4  }, {0,0,1,2,3  }, x, v);
  case 12: return evalon<12,6>(xnode, {4,4,4,4,4,4}, {0,0,1,2,3,4}, x, v);
  case 13: return evalon<13,6>(xnode, {4,4,4,4,4,4}, {0,0,1,2,3,4}, x, v);
  }
  return false;
}

bool UniformOffsetNodal
::eval_derivative (const Int& np, const Real& x, Real* const v) const {
  assert(0);
  return false;  
}

bool FreeNodal::get_x (const Int& np, const Real*& x) const {
  switch (np) {
  case 4: x = x_np4; break;
  case 5: x = x_np5; break;
  case 6: x = x_np6; break;
  case 7: x = x_np7; break;
  case 8: x = x_np8; break;
  case 10: x = x_np10; break;
  default: x = nullptr; return false;
  }
  return true;
}

bool FreeNodal::get_w (const Int& np, const Real*& w) const {
  switch (np) {
  case 4: w = w_np4; break;
  case 5: w = w_np5; break;
  case 6: w = w_np6; break;
  case 7: w = w_np7; break;
  case 8: w = w_np8; break;
  case 10: w = w_np10; break;
  default: w = nullptr; return false;
  }
  return true;
}

bool FreeNodal::eval (const Int& np, const Real& x, Real* const v) const {
  const Real* xnodes;
  get_x(np, xnodes);
  switch (np) {
  case 4: {
    const Int subnp[] = {3,4};
    const Int all[] = {0,1,2,3};
    const Int  n0[] = {0,1,2  };
    const Int* nodes[] = {n0, all};
    islet::eval(4, xnodes, subnp, nodes, x, v);
  } break;
  case 5: {
    const Int subnp[] = {4,4};
    const Int n01[] = {0,1,2,3};
    const Int* nodes[] = {n01, n01};
    islet::eval(5, xnodes, subnp, nodes, x, v);
  } break;
  case 6: {
    const Int subnp[] = {5,5,6};
    const Int all[] = {0,1,2,3,4,5};
    const Int n01[] = {0,1,2,3,4, };
    const Int* nodes[] = {n01, n01, all};
    islet::eval(6, xnodes, subnp, nodes, x, v);
  } break;
  case 7: {
    const Int subnp[] = {6,7,6};
    const Int all[] = {0,1,2,3,4,5,6};
    const Int  n0[] = {0,1,2,  4,5,6};
    const Int  n2[] = {0,1,2,3,4,5  };
    const Int* nodes[] = {n0, all, n2};
    islet::eval(7, xnodes, subnp, nodes, x, v);
  } break;
  case 8: {
    const Int subnp[] = {7,8,7,8};
    const Int all[] = {0,1,2,3,4,5,6,7};
    const Int  n0[] = {0,1,2,3,4,5,  7};
    const Int  n2[] = {0,1,2,3,4,5,6  };
    const Int* nodes[] = {n0, all, n2, all};
    islet::eval(8, xnodes, subnp, nodes, x, v);
  } break;
  case 10: {
    const Int subnp[] = {9,9,10,9,10};
    const Int all[] = {0,1,2,3,4,5,6,7,8,9};
    const Int  n0[] = {0,1,2,3,4,5,6,7,8  };
    const Int  n1[] = {0,1,2,3,4,5,  7,8,9};
    const Int  n3[] = {0,1,2,3,4,5,6,7,8  };
    const Int* nodes[] = {n0, n1, all, n3, all};
    islet::eval(10, xnodes, subnp, nodes, x, v);
  } break;
  default: return false;
  }
  return true;
}

bool FreeNodal::eval_derivative (const Int& np, const Real& x, Real* const v) const {
  // For p-refinement, I believe we do not need this function. Assert that. If
  // we do, I need to think about how I want to handle it, since these basis
  // functions have discontinuous first derivs at the nodes.
  assert(0);
  return false;
}

} // namespace islet
} // namespace slmm
