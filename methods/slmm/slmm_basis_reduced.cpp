#include "slmm_basis_reduced.hpp"

namespace slmm {
using namespace islet;

bool UniformNodeReduced::get_x (const Int& np, const Real*& x) const {
  if ((np > 13 && np < GLL::np_max) || np > GLL::np_max) return false;
  static Real xs[((GLL::np_max+1)*GLL::np_max)/2] = {0};
  x = xs + ((np-1)*np)/2;
  if (x[0] == 0) {
    Real* const xsub = const_cast<Real*>(x);
    xsub[0] = -1; xsub[np-1] = 1;
    for (Int j = 1; j < np-1; ++j) xsub[j] = 2*(Real(j)/(np-1))-1;
  }
  return true;
}

bool UniformNodeReduced::get_w (const Int& np, const Real*& w) const {
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

bool UniformNodeReduced
::eval_derivative (const Int& np, const Real& x, Real* const v) const {
  assert(0);
  return false;
}

Int UniformNodeReduced::max_degree (const Int& np) const { return 1; }

bool UniformNodeReduced::eval (const Int& np, const Real& x, Real* const v) const {
  const Real* xnode;
  get_x(np, xnode);
  switch (np) {
  case  2: return evalon< 2,0>(xnode, {           }, {           }, x, v);
  case  3: return evalon< 3,1>(xnode, {2          }, {0          }, x, v);
  case  4: return evalon< 4,2>(xnode, {2,2        }, {0,1        }, x, v);
  case  5: return evalon< 5,2>(xnode, {2,2        }, {0,1        }, x, v);
  case  6: return evalon< 6,3>(xnode, {2,2,2      }, {0,1,2      }, x, v);
  case  7: return evalon< 7,3>(xnode, {2,2,2      }, {0,1,2      }, x, v);
  case  8: return evalon< 8,4>(xnode, {2,2,2,2    }, {0,1,2,3    }, x, v);
  case  9: return evalon< 9,4>(xnode, {2,2,2,2    }, {0,1,2,3    }, x, v);
  case 10: return evalon<10,5>(xnode, {2,2,2,2,2  }, {0,1,2,3,4  }, x, v);
  case 11: return evalon<11,5>(xnode, {2,2,2,2,2  }, {0,1,2,3,4  }, x, v);
  case 12: return evalon<12,6>(xnode, {2,2,2,2,2,2}, {0,1,2,3,4,5}, x, v);
  case 13: return evalon<13,6>(xnode, {2,2,2,2,2,2}, {0,1,2,3,4,5}, x, v);
  case 16: return evalon<16,8>(xnode, {2,2,2,2,2,2,2,2}, {0,1,2,3,4,5,6,7}, x, v);
  default: return false;
  }
  return true;
}

} // namespace slmm
