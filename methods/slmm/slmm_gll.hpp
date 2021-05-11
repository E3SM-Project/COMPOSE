#ifndef INCLUDE_SLMM_GLL_HPP
#define INCLUDE_SLMM_GLL_HPP

#include "slmm_basis.hpp"

namespace slmm {

class GLL : public Basis {
  const Real oo3 = 1.0/3.0;
  const Real to3 = 2.0/3.0;
  const Real oo6 = 1.0/6.0;
  const Real oo8 = 1.0/8.0;
  const Real sqrt5 = std::sqrt(5.0);
  const Real oosqrt5 = 1.0/sqrt5;
  const Real np2_coord[2] = {-1.0, 1.0};
  const Real np2_wt[2]    = {1.0, 1.0};
  const Real np3_coord[3] = {-1.0, 0.0, 1.0};
  const Real np3_wt[3]    = {oo3, 2.0 - to3, oo3};
  const Real np4_coord[4] = {-1.0, -oosqrt5, oosqrt5, 1.0};
  const Real np4_wt[4]    = {oo6, 1.0 - oo6, 1.0 - oo6, oo6};
#if SLMM_NP_MAX >= 5
  const Real sqrt3o7 = std::sqrt(3.0/7.0);
  const Real np5_coord[5] = {-1.0, -sqrt3o7, 0, sqrt3o7, 1.0};
  const Real np5_wt[5]    = {0.1, 49.0/90.0, 32.0/45.0, 49.0/90.0, 0.1};
#endif
#if SLMM_NP_MAX >= 6
  const Real np6_coord[6] = {
# define np6a std::sqrt(1.0/3.0 + 2.0*std::sqrt(7.0)/21.0)
# define np6b std::sqrt(1.0/3.0 - 2.0*std::sqrt(7.0)/21.0)
    -1, -np6a, -np6b, np6b, np6a, 1
# undef np6b
# undef np6a
  };
  const Real np6_wt[6] = {
# define v0 1.0/15.0
# define v1 (14 - std::sqrt(7.0))/30.0
# define v2 (14 + std::sqrt(7.0))/30.0
    v0, v1, v2, v2, v1, v0
# undef v2
# undef v1
# undef v0
  };
#endif
#if SLMM_NP_MAX >= 7
  const Real np7_coord[7] = {
# define np7a std::sqrt((5.0 + 2.0*std::sqrt(5.0/3.0))/11.0)
# define np7b std::sqrt((5.0 - 2.0*std::sqrt(5.0/3.0))/11.0)
    -1, -np7a, -np7b, 0, np7b, np7a, 1,
# undef np7b
# undef np7a
  };
  const Real np7_wt[7] = {
# define v0 1.0/21.0
# define v1 (124 - 7*std::sqrt(15.0))/350.0
# define v2 (124 + 7*std::sqrt(15.0))/350.0
    v0, v1, v2, 256.0/525.0, v2, v1, v0,
# undef v0
# undef v1
# undef v2
  };
#endif
#if SLMM_NP_MAX >= 8
  const Real np8_coord[8] = {
    -1, -0.8717401485096066153, -0.59170018143314230214, -0.20929921790247886877,
    0.20929921790247886877, 0.59170018143314230214, 0.87174014850960661534, 1
  };
  const Real np8_wt[8] = {
    0.03571428571428571429, 0.21070422714350603938, 0.34112269248350436476,
    0.41245879465870388157, 0.41245879465870388157, 0.34112269248350436476,
    0.21070422714350603938, 0.03571428571428571429
  };
#endif
#if SLMM_NP_MAX >= 9
const Real np9_coord[9] = {
    -1, -0.89975799541146015731, -0.67718627951073775345, -0.36311746382617815871,
    0, 0.36311746382617815871, 0.67718627951073775345, 0.89975799541146015731, 1
  };
  const Real np9_wt[9] = {
    0.02777777777777777778, 0.16549536156080552505, 0.27453871250016173528,
    0.34642851097304634512, 0.37151927437641723356, 0.34642851097304634512,
    0.27453871250016173528, 0.16549536156080552505, 0.02777777777777777778
  };
#endif
#if SLMM_NP_MAX >= 10
const Real np10_coord[10] = {
  -1, -0.91953390816645881383, -0.73877386510550507500, -0.47792494981044449566,
  -0.16527895766638702463, 0.16527895766638702463, 0.47792494981044449566,
  0.73877386510550507500, 0.91953390816645881383, 1
};
const Real np10_wt[10] = {
  0.02222222222222222222, 0.13330599085107011113, 0.22488934206312645212,
  0.29204268367968375788, 0.32753976118389745666, 0.32753976118389745666,
  0.29204268367968375788, 0.22488934206312645212, 0.13330599085107011113,
  0.02222222222222222222
};
#endif
#if SLMM_NP_MAX >= 11
const Real np11_coord[11] = {
  -1, -0.93400143040805913433, -0.78448347366314441862, -0.56523532699620500647,
  -0.29575813558693939143, 0, 0.29575813558693939143, 0.56523532699620500647,
    0.78448347366314441862, 0.93400143040805913433, 1
  };
  const Real np11_wt[11] = {
    0.01818181818181818182, 0.10961227326699486446, 0.18716988178030520411,
    0.24804810426402831404, 0.28687912477900808868, 0.30021759545569069379,
    0.28687912477900808868, 0.24804810426402831404, 0.18716988178030520411,
    0.10961227326699486446, 0.01818181818181818182
  };
#endif
#if SLMM_NP_MAX >= 12
  const Real np12_coord[12] = {
    -1, -0.94489927222288222341, -0.81927932164400667835, -0.63287615303186067766,
    -0.39953094096534893226, -0.13655293285492755486, 0.13655293285492755486,
    0.39953094096534893226, 0.63287615303186067766, 0.81927932164400667835,
    0.94489927222288222341, 1
  };
  const Real np12_wt[12] = {
    0.01515151515151515152, 0.09168451741319613067, 0.15797470556437011517,
    0.21250841776102114536, 0.25127560319920128029, 0.27140524091069617700,
    0.27140524091069617700, 0.25127560319920128029, 0.21250841776102114536,
    0.15797470556437011517, 0.09168451741319613067, 0.01515151515151515152
  };
#endif
#if SLMM_NP_MAX >= 13
  const Real np13_coord[13] = {
    -1, -0.95330984664216391190, -0.84634756465187231687, -0.68618846908175742607,
    -0.48290982109133620175, -0.24928693010623999257, 0, 0.24928693010623999257,
    0.48290982109133620175, 0.68618846908175742607, 0.84634756465187231687,
    0.95330984664216391190, 1
  };
  const Real np13_wt[13] = {
    0.01282051282051282051, 0.07780168674681892779, 0.13498192668960834912,
    0.18364686520355009201, 0.22076779356611008609, 0.24401579030667635646,
    0.25193084933344673604, 0.24401579030667635646, 0.22076779356611008609,
    0.18364686520355009201, 0.13498192668960834912, 0.07780168674681892779,
    0.01282051282051282051
  };
#endif
#if SLMM_NP_MAX >= 16
  const Real np16_coord[16] = {
    -1, -0.96956804627021793295, -0.89920053309347209299, -0.79200829186181506393,
    -0.65238870288249308947, -0.48605942188713761178, -0.29983046890076320810,
    -0.10132627352194944784, 0.10132627352194944784, 0.29983046890076320810,
    0.48605942188713761178, 0.65238870288249308947, 0.79200829186181506393,
    0.89920053309347209299, 0.96956804627021793295, 1
  };
  const Real np16_wt[16] = {
    0.00833333333333333333, 0.05085036100591990540, 0.08939369732593080099,
    0.12425538213251409835, 0.15402698080716428081, 0.17749191339170412530,
    0.19369002382520358432, 0.20195830817822987149, 0.20195830817822987149,
    0.19369002382520358432, 0.17749191339170412530, 0.15402698080716428081,
    0.12425538213251409835, 0.08939369732593080099, 0.05085036100591990540,
    0.00833333333333333333
  };
#endif

public:

  const char* name () const override { return "GLL"; }

  KOKKOS_INLINE_FUNCTION GLL () {}

  KOKKOS_INLINE_FUNCTION
  bool get_coef (const Int& np, const Real*& coord, const Real*& wt) const {
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
#if SLMM_NP_MAX >= 5
    case 5:
      coord = np5_coord;
      wt = np5_wt;
      break;
#endif
#if SLMM_NP_MAX >= 6
    case 6:
      coord = np6_coord;
      wt = np6_wt;
      break;
#endif
#if SLMM_NP_MAX >= 7
    case 7:
      coord = np7_coord;
      wt = np7_wt;
      break;
#endif
#if SLMM_NP_MAX >= 8
    case 8:
      coord = np8_coord;
      wt = np8_wt;
      break;
#endif
#if SLMM_NP_MAX >= 9
    case 9:
      coord = np9_coord;
      wt = np9_wt;
      break;
#endif
#if SLMM_NP_MAX >= 10
    case 10:
      coord = np10_coord;
      wt = np10_wt;
      break;
#endif
#if SLMM_NP_MAX >= 11
    case 11:
      coord = np11_coord;
      wt = np11_wt;
      break;
#endif
#if SLMM_NP_MAX >= 12
    case 12:
      coord = np12_coord;
      wt = np12_wt;
      break;
#endif
#if SLMM_NP_MAX >= 13
    case 13:
      coord = np13_coord;
      wt = np13_wt;
      break;
#endif
#if SLMM_NP_MAX >= 16
    case 16:
      coord = np16_coord;
      wt = np16_wt;
      break;
#endif
    default:
      coord = nullptr;
      wt = nullptr;
      return false;
    }
    return true;
  }

  bool get_x (const Int& np, const Real*& coord) const override {
    const Real* wt;
    return get_coef(np, coord, wt);
  }

  bool get_w (const Int& np, const Real*& wt) const override {
    const Real* coord;
    return get_coef(np, coord, wt);
  }

  // x in [-1, 1].
  KOKKOS_INLINE_FUNCTION
  bool eval (const Int& np, const Real& x, Real* const ge) const override {
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
      const Real x2 = x*x;
      ge[0] = (1.0 - x)*(5.0*x2 - 1.0)*oo8;
      ge[1] = -sqrt5*oo8*(sqrt5 - 5.0*x)*(x2 - 1.0);
      ge[2] = -sqrt5*oo8*(sqrt5 + 5.0*x)*(x2 - 1.0);
      ge[3] = (1.0 + x)*(5.0*x2 - 1.0)*oo8;
    } break;
    default: {
#if SLMM_NP_MAX > 4
      if (np > SLMM_NP_MAX) ko::abort("GLL::eval: order not supported.");
      const Real* coord, * wt;
      get_coef(np, coord, wt);
      eval_lagrange_poly<Real>(np, coord, x, ge);
      break;
#endif
      for (int i = 0; i < np; ++i) ge[i] = 0;
      return false;
    }
    }
    return true;
  }

  KOKKOS_INLINE_FUNCTION
  bool eval_derivative (const Int& np, const Real& x, Real* const ge) const override {
    switch (np) {
    case 2: {
      ge[0] = -0.5;
      ge[1] =  0.5;
    } break;
    case 3: {
      const Real x2p = 2*x;
      ge[0] =  0.5*(x2p - 1);
      ge[1] = -x2p;
      ge[2] =  0.5*(x2p + 1);
    } break;
    case 4: {
      ge[0] =  oo8*((10 - 15*x)*x + 1);
      ge[1] = -sqrt5*oo8*((2*sqrt5 - 15*x)*x + 5);
      ge[2] = -sqrt5*oo8*((2*sqrt5 + 15*x)*x - 5);
      ge[3] =  oo8*((10 + 15*x)*x - 1);
    } break;
    default:
#if SLMM_NP_MAX > 4
      if (np > SLMM_NP_MAX) ko::abort("GLL::eval_derivative: order not supported.");
      const Real* coord, * wt;
      get_coef(np, coord, wt);
      eval_lagrange_poly_derivative<Real>(np, coord, x, ge);
      break;
#endif
      for (int i = 0; i < np; ++i) ge[i] = 0;
      return false;
    }
    return true;
  }
};

} // namespace slmm

#endif
