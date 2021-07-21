#include "islet_npx.hpp"
#include "islet_util.hpp"

void op_eval (const InterpMethod& im, const Real a_src, Real* v) {
  switch (im.type) {
  case InterpMethod::npx: npx<Real>::eval(im.np, a_src, v); break;
  case InterpMethod::npxstab: npxstab<Real>::eval(im.np, a_src, v); break;
  case InterpMethod::user: im.uim->eval(a_src, v); break;
  default:
    throw_if(true, "op_eval: invalid im.type " << im.type);
  }
}

template <typename Scalar>
void npxstab<Scalar>::eval (const Int& np, const Scalar& x, Scalar* const v) {
  switch (np) {                                                   // order of accuracy
  case  2: return eval< 2,0>({                      }, {                }, x, v); // 1
  case  3: return eval< 3,0>({                      }, {                }, x, v); // 2
  case  4: return eval< 4,1>({ 3                    }, {0,              }, x, v); // 2
  case  5: return eval< 5,2>({ 3,  4                }, {0, 0            }, x, v); // 2
  case  6: return eval< 6,2>({ 6,  5                }, {0, 0            }, x, v); // 4
  case  7: return eval< 7,3>({ 5,  5,  6            }, {0, 0, 0         }, x, v); // 4
  case  8: return eval< 8,4>({ 6,  6,  7,  6        }, {0, 0, 0, 1      }, x, v); // 5
  case  9: return eval< 9,4>({ 7,  8,  7,  7        }, {0, 0, 0, 1      }, x, v); // 6
  case 10: return eval<10,5>({ 7,  7,  7,  8,  8    }, {0, 0, 0, 0, 1   }, x, v); // 6
  case 11: return eval<11,5>({ 8,  9,  8,  9,  8    }, {0, 0, 0, 0, 1   }, x, v); // 7
  case 12: return eval<12,6>({ 9,  9, 10, 10,  9, 10}, {0, 0, 0, 0, 1, 1}, x, v); // 8
  case 13: return eval<13,6>({10, 10, 10, 10, 11, 10}, {0, 0, 0, 0, 0, 1}, x, v); // 9
  case 16: return eval<16,8>({12, 13, 13, 13, 13, 14, 13, 12},
                              { 0,  0,  0,  0,  0,  0,  1,  2}, x, v); // 11
  default: throw_if(true, "Only 2 <= np <= 13, np = 16 are supported.");
  }
}

template class npxstab<Real>;
