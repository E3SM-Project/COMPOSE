#ifndef INCLUDE_MPAS_ORDER_INVARIANT_ACCUM_INTREPS_HPP
#define INCLUDE_MPAS_ORDER_INVARIANT_ACCUM_INTREPS_HPP

#include <type_traits>

#include "mpas_util.hpp"

namespace order_invariant_accum {

// Helper for atomic (real, integer) pairs.
template <typename RealT, typename IntRep, bool atomic_accum>
struct RealIntRep {
  typedef RealT Real;

  RealIntRep () : v(0) {}
  RealIntRep (const Real a) { v = mult*a; }

  Real real () const { return Real(v)/mult; }

  RealIntRep& operator+= (const RealIntRep& a) {
    if (atomic_accum) ompatomic v += a.v;
    else v += a.v;
    return *this;
  }

  bool is_zero () const { return v == 0; }

private:
  static constexpr IntRep mult = std::numeric_limits<IntRep>::max();
  IntRep v;
};

template <typename RealT, typename Int>
struct RealInt2Rep {
  typedef RealT Real;

  RealInt2Rep () : v({0,0}) {}
  RealInt2Rep (const Real a) { d2i(a, v); }

  Real real () const { return i2d(v); }

  RealInt2Rep& operator+= (const RealInt2Rep& a) {
    sum_into(a.v, v);
    return *this;
  }

  bool is_zero () const { return v[0] == 0 && v[1] == 0; }

private:
  typedef typename std::make_unsigned<Int>::type UInt;
  static constexpr int nbitm1 = 8*sizeof(Int) - 1; // number of bits minus one
  static constexpr UInt mult = std::numeric_limits< Int>::max(); // 01...1
  typedef std::array<UInt,2> A;

  template <typename T> static inline T intbits (const T& v) { return  mult & v; }
  template <typename T> static inline T signbit (const T& v) { return ~mult & v; }

  // This impl uses two's complement in a list (hard-set to two) of unsigned
  // integers. The 0'th integer has the sign bit as its leftmost bit. The other
  // integers (just one in this impl) use all but the leftmost bit to store
  // data. In sum_into, the leftmost bit can then be used as the carry bit.

  static inline void d2i (const Real a, A& is) {
    assert(std::abs(a) <= 1);
    const Real b = std::abs(mult*a);
    is[0] = std::min(mult, UInt(b));
    is[1] = mult*(b - is[0]);
    if (a < 0) {
      const UInt n1 = intbits(~is[1]) + 1; // 2's compl
      is[1] = intbits(n1);
      is[0] = ~is[0] + (n1 >> nbitm1);     // 2's compl and carry
      assert(signbit(is[0]) || (is[0] == 0 && is[1] == 0));
    }
    assert( ! signbit(is[1]));
  }

  static inline Real i2d (const A& is) {
    if (signbit(is[0])) {
      const UInt is0 = intbits(~is[0]) + 1;
      const UInt is1 = intbits(~is[1]) + (is0 >> nbitm1);
      const auto v = -(Real(is0  ) + Real(is1  )/mult)/mult;
      assert(v <  0);
      return v;
    } else {
      const auto v =  (Real(is[0]) + Real(is[1])/mult)/mult;
      assert(v >= 0);
      return v;
    }
  }

  static inline void sum_into (const A& s, A& d) {
    const UInt sum = s[1] + d[1];
    d[1]  = intbits(sum);
    d[0] += s[0] + (sum >> nbitm1); // second term is the carry bit
  }

  friend bool unittest_RealInt2Rep();

private:
  A v;
};

template <typename Real, bool extended_precision, bool atomic_accum>
struct IntRep;

typedef RealTypeTraits<double>::SameSizeInt DInt;
typedef RealTypeTraits<float>::SameSizeInt FInt;

template <> struct IntRep<double,false,false>
  : public RealIntRep<double,DInt,false> {
  typedef RealIntRep<double,DInt,false> Super;
  typedef DInt Int;
  IntRep () : Super() {}
  IntRep (const double a) : Super(a) {}
};

template <> struct IntRep<double,false,true>
  : public RealIntRep<double,DInt,true> {
  typedef RealIntRep<double,DInt,true> Super;
  typedef DInt Int;
  IntRep () : Super() {}
  IntRep (const double a) : Super(a) {}
};

template <> struct IntRep<double,true,false>
  : public RealInt2Rep<double,DInt> {
  typedef RealInt2Rep<double,DInt> Super;
  typedef DInt Int;
  IntRep () : Super() {}
  IntRep (const double a) : Super(a) {}
};

template <> struct IntRep<float,true,false>
  : public RealIntRep<float,DInt,false> {
  typedef RealIntRep<float,DInt,false> Super;
  typedef DInt Int;
  IntRep () : Super() {}
  IntRep (const float a) : Super(a) {}
};

template <> struct IntRep<float,false,true>
  : public RealIntRep<float,FInt,true> {
  typedef RealIntRep<float,FInt,true> Super;
  typedef DInt Int;
  IntRep () : Super() {}
  IntRep (const float a) : Super(a) {}
};

template <> struct IntRep<float,true,true>
  : public RealIntRep<float,DInt,true> {
  typedef RealIntRep<float,DInt,true> Super;
  typedef DInt Int;
  IntRep () : Super() {}
  IntRep (const float a) : Super(a) {}
};

bool unittest_RealInt2Rep();

} // namespace order_invariant_accum

#endif
