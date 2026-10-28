#include "mpas_order_invariant_accum.hpp"
#include "mpas_util.hpp"

namespace order_invariant_accum {

// Test that this gives nearly the same results as IntRep<float, true,
// false>. Can't be exactly the same b/c it has one fewer bit than the other.
// There is no reason to use these in practice.
template <> struct IntRep<float,false,false>
  : public RealIntRep<float,int32_t,false> {
  typedef RealIntRep<float,int32_t,false> Super;
  typedef int32_t Int;
  IntRep () : Super() {}
  IntRep (const float a) : Super(a) {}
};

struct FloatInt2Rep : public RealInt2Rep<float,int32_t> {
  typedef RealInt2Rep<float,int32_t> Super;
  FloatInt2Rep (const float a) : Super(a) {}
};

template <typename Real, bool extended_precision, bool atomic_accum,
          bool use_FloatInt2Rep=false>
struct TestAccumType {
  typedef typename AccumType<Real,extended_precision,atomic_accum>::type type;
};

template <typename Real>
struct TestAccumType<Real, true, false, true> {
  typedef Accumulator<FloatInt2Rep> type;
};

static uint32_t intbits (const uint32_t& v) { return ((1u << 31)-1) & v; }
static uint32_t signbit (const uint32_t& v) { return  (1u << 31)    & v; }

static std::array<uint32_t,2> I2ii (const int64_t a) {
  if (a < 0) {
    const uint64_t aus = -a;
    const uint32_t i1 = intbits(~intbits(aus)) + 1;
    return std::array<uint32_t,2>(
      {(1u << 31) | (intbits(~uint32_t(aus >> 31)) + (i1 >> 31)),
       intbits(i1)});
  } else {
    const uint64_t aus = a;
    return std::array<uint32_t,2>({uint32_t(aus >> 31),
                                   uint32_t(intbits(aus))});
  }
}

static int64_t ii2I (const std::array<uint32_t,2>& a2) {
  if (signbit(a2[0])) {
    const uint32_t i1 = intbits(~a2[1]) + 1;
    return -int64_t((uint64_t(intbits(~a2[0]) + (i1 >> 31)) << 31) |
                    uint64_t(intbits(i1)));
  } else {
    return int64_t((uint64_t(a2[0]) << 31) | uint64_t(a2[1]));
  }
}

static std::string to_string (const std::array<uint32_t,2>& a2) {
  std::stringstream ss;
  ss << "(" << a2[0] << ", " << a2[1] << ")";
  return ss.str();
}

bool unittest_RealInt2Rep () {
  typedef RealInt2Rep<float,int32_t> S;
  typedef RealInt2Rep<double,int64_t> D;
  typedef int64_t lint;
  bool ok = true;
  // Test integer arithmetic.
  const auto f = [&] (const lint a, const lint b) {
    const lint c = a + b;
    S::A a2 = I2ii(a), b2 = I2ii(b);
    if (!(ii2I(a2) == a && ii2I(b2) == b)) ok = false;
    S::A c2 = a2;
    assert(ii2I(c2) == a);
    S::sum_into(b2, c2);
    if (ii2I(c2) != c) {
      const auto c_check = I2ii(c);
      printf("S::sum_into FAILed:");
      printf("  a: %ld %s b: %ld %s\n",
             a, to_string(a2).c_str(), b, to_string(b2).c_str());
      printf("  c: %ld (%d, %d) %ld (%d, %d)\n",
             c, c_check[0], c_check[1], ii2I(c2), c2[0], c2[1]);
      ok = false;
    }
  };
  {
    const lint
      a = (lint(S::mult) << 5) + (12345 << 4),
      b = (lint(S::mult) << 12) + (983987654321 << 3);
    f(b,a); f(b,-a); f(-b,a); f(-b,-a);
    f(a,b); f(a,-b); f(-a,b); f(-a,-b);
  }
  {
    const lint
      a = (lint(S::mult) << 5) + (983987654321 << 3),
      b = (lint(S::mult) << 12) + (12345 << 4);
    f(b,a); f(b,-a); f(-b,a); f(-b,-a);
    f(a,b); f(a,-b); f(-a,b); f(-a,-b);
  }
  // Test real <-> integer rep.
  {
    static const float  stol = std::numeric_limits<float >::epsilon();
    static const double dtol = std::numeric_limits<double>::epsilon();
    const double c = 0.1*M_PI;
    if (!(std::abs(S(c).real() - c)  <= c*stol)) ok = false;
    if (!(std::abs(D(c).real() - c)  <= c*dtol)) ok = false;
    if (!(std::abs(S(-c).real() + c) <= c*stol)) ok = false;
    if (!(std::abs(D(-c).real() + c) <= c*dtol)) ok = false;
  }
  return ok;
}

struct TestOutput {
  double sum, re;
  TestOutput () : sum(0), re(-1) {}
};

template <typename Real, bool extended_precision, bool atomic_accum,
          bool use_FloatInt2Rep=false>
static bool test (const int n, const bool general,
                  const int seed = -1, TestOutput* out = nullptr) {
  const bool verbose = false;
  typedef typename order_invariant_accum
    ::TestAccumType<Real, extended_precision, atomic_accum, use_FloatInt2Rep>::type
    Accum;
  if (verbose)
    printf("%s%s%s%s\n",
           sizeof(Real) == sizeof(float) ? "float" : "double",
           extended_precision ? " extended_precision" : "",
           atomic_accum ? " atomic_accum" : "",
           use_FloatInt2Rep ? " FloatInt2Rep" : "");
  bool ok = true;
  // Fill terms to sum.
  std::vector<Real> v(n);
  const Real os = general ? -0.5 : 0;
  if (seed >= 0) srand(seed);
  for (int i = 0; i < n; ++i)
    v[i] = (urand() + os)*std::pow(urand(), 3);
  for (int i = 0; i < n; ++i) v[i] *= -1;
  // Compute max term magnitude.
  Real maxmag = 0;
  for (int i = 0; i < n; ++i)
    maxmag = std::max(maxmag, std::abs(v[i]));
  Real sum0 = 0;
  int nout = 0;
  const int nthr = mpas_omp_get_max_threads();
  Ar1 wrk(nthr*Accum::nreal);
  for (int trial = 0; trial < std::min(10*n, 10000); ++trial) {
    // Determine a new permutation.
    for (int i = 0; i < std::min(1000, n); ++i) {
      const int i0 = int(urand()*n);
      const int i1 = int(urand()*n);
      std::swap(v[i0], v[i1]);
    }
    // Sum v, testing both accumulation ops and both 'sum' routines.
    auto s = Accum::sum(wrk, n/2, v.data(), n*maxmag);
    Accum::sum(wrk, v.data() + n/2, n - n/2, &s, 1, true);
    // Check result.
    const Real sum = s();
    if (trial == 0) {
      sum0 = sum;
      Real rsum = 0;
      for (int i = 0; i < n; ++i)
        rsum += v[i];
      const auto re = std::abs((sum - rsum)/rsum);
      if (verbose)
        printf("maxmag %23.16e rsum %23.16e diff %23.16e re %9.2e\n",
               maxmag, rsum, sum - rsum, re);
      if (out) {
        out->sum = sum;
        out->re = re;
      }
      if ( ! general && re > n*std::numeric_limits<Real>::epsilon())
        ok = false;
    } else {
      if (sum != sum0) {
        ok = false;
        if (nout < 10) {
          printf("trial %4d sum0 %23.16e sum %10.3e diff %10.3e re %9.2e\n",
                 trial, sum0, sum, sum - sum0, (sum - sum0)/sum0);
          ++nout;
        }
      }
    }
  }
  if ( ! ok) printf("Accumulator test FAILed\n");
  return ok;
}

bool unittest_Accumulator () {
  bool ok = true;
  
  if ( ! order_invariant_accum::unittest_RealInt2Rep()) {
    ok = false;
    printf("order_invariant_accum::unittest_RealInt2Rep FAILed\n");
  }

  const int n = 999;

  const int seed = 42;
  TestOutput o1t, o2t, o1f;
  test<float,true ,false     >(n, true, seed, &o1t);
  test<float,true ,false,true>(n, true, seed, &o2t);
  test<float,false,false     >(n, true, seed, &o1f);
  if (o2t.re > o1f.re) {
    printf("Int2 vs int64_t FAILed with larger re\n");
    ok = false;
  }
  if (o1f.re/o2t.re < 0.25*o1f.re/o1t.re) {
    printf("Int2 vs int64_t FAILed with too small of a re gain:"
           " %f vs %f\n", o1f.re/o2t.re, o1f.re/o1t.re);
    ok = false;
  }

  if ( ! test<float ,false,false>(n, false)) ok = false;
  if ( ! test<float ,true ,false>(n, false)) ok = false;
  if ( ! test<float ,true ,true >(n, true )) ok = false;
  if ( ! test<float ,true ,true >(n, false)) ok = false;
  if ( ! test<double,false,true >(n, true )) ok = false;
  if ( ! test<double,false,false>(n, false)) ok = false;
  if ( ! test<double,true ,false>(n, true )) ok = false;
  if ( ! test<double,true ,false>(n, false)) ok = false;

  return ok;
}

} // namespace order_invariant_accum
