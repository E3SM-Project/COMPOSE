#ifndef INCLUDE_SLMM_GLL_HPP
#define INCLUDE_SLMM_GLL_HPP

#include "slmm_defs.hpp"

namespace slmm {

class GLL {
  const Real oo3 = 1.0/3.0;
  const Real to3 = 2.0/3.0;
  const Real sqrt5 = std::sqrt(5.0);
  const Real oo6 = 1.0/6.0;
  const Real np2_coord[2] = {-1.0, 1.0};
  const Real np2_wt[2]    = {1.0, 1.0};
  const Real np3_coord[3] = {-1.0, 0.0, 1.0};
  const Real np3_wt[3]    = {oo3, 2.0 - to3, oo3};
  const Real np4_coord[4] = {-1.0, -1.0/sqrt5, 1.0/sqrt5, 1.0};
  const Real np4_wt[4]    = {oo6, 1.0 - oo6, 1.0 - oo6, oo6};

public:
  enum { max_np = 4 };

  KOKKOS_INLINE_FUNCTION GLL () {}

  KOKKOS_INLINE_FUNCTION
  void get_coef (const int np, const Real*& coord, const Real*& wt) {
    switch (np) {
    case 2:
      coord = np2_coord;
      wt = np2_wt;
      break;
    case 3:
      coord = np3_coord;
      wt = np3_wt;
      break;
    case 4:
      coord = np4_coord;
      wt = np4_wt;
      break;
    default:
      ko::abort("GLL::get_coef: order not supported.");
    }
  }

  // x in [-1, 1].
  KOKKOS_INLINE_FUNCTION
  void eval (const int np, const Real& x, Real* const ge) const {
    switch (np) {
    case 2: {
      ge[0] = 0.5*(1.0 - x);
      ge[1] = 0.5*(1.0 + x);
    } break;
    case 3: {
      const Real x2 = x*x;
      ge[0] = 0.5*(x2 - x);
      ge[1] = 1.0 - x2;
      ge[2] = 0.5*(x2 + x);
    } break;
    case 4: {
      const Real oo8 = 1.0/8.0;
      const Real x2 = x*x;
      ge[0] = (1.0 - x)*(5.0*x2 - 1.0)*oo8;
      ge[1] = -sqrt5*oo8*(sqrt5 - 5.0*x)*(x2 - 1.0);
      ge[2] = -sqrt5*oo8*(sqrt5 + 5.0*x)*(x2 - 1.0);
      ge[3] = (1.0 + x)*(5.0*x2 - 1.0)*oo8;
    } break;
    default:
      ko::abort("GLL::eval: order not supported.");
    }
  }
};

} // namespace slmm

#endif
