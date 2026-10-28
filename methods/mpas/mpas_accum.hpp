#ifndef INCLUDE_MPAS_ACCUM_HPP
#define INCLUDE_MPAS_ACCUM_HPP

#include "mpas_types.hpp"
#include "mpas_order_invariant_accum.hpp"

namespace accum {

static const Int threaded_bfb_nslot = 512;

void reduce_in_place(Real* a, const Int nscalar, const Int n);

// Sum-reduce fn(i, accum) for i in [0,n). accum[k] is 0-inited for k in
// [0,nscalar). In fn, update accum with =, not +=.
template <int max_nscalar, typename Fn>
void accumulate_threaded_bfb (Fn fn, const Int n, Real* accum,
                              const Int nscalar = max_nscalar) {
  static const Int nslot = threaded_bfb_nslot;
  Real slot[nslot*max_nscalar];
  const Int nperslot = (n + nslot - 1)/nslot;
  ompparfor for (Int si = 0; si < nslot; ++si) {
    Real* slotp = slot + si*nscalar;
    for (Int j = 0; j < nscalar; ++j) slotp[j] = 0;
    const Int ib = si*nperslot;
    const Int ie = std::min(ib + nperslot, n);
    for (Int i = ib; i < ie; ++i) {
      Real tmp[max_nscalar];
      fn(i, tmp);
      for (Int j = 0; j < nscalar; ++j) slotp[j] += tmp[j];
    }
  }
  reduce_in_place(slot, nscalar, threaded_bfb_nslot);
  for (Int j = 0; j < nscalar; ++j) accum[j] = slot[j];
}

Real accumulate_threaded_bfb(Real* const a, const Int n,
                             Real slot[threaded_bfb_nslot],
                             const bool in_threaded_region = true);

} // namespace accum

/* Algorithm from https://en.wikipedia.org/wiki/Kahan_summation_algorithm
     This is a compensated accumulator. If every element in the list of n
   elements to sum has the same sign, then the condition number of the sum is 1
   and so the error bound on the output is O(2 eps + n eps^2). Thus, if n <
   ~1/e, then the output is independent of summation order, up to 2 eps.
 */
struct KahanBabushkaNeumaierSummer {
  typedef Real scalar_type;

  KahanBabushkaNeumaierSummer () : sum(0), c(0) {}
  KahanBabushkaNeumaierSummer (const Real accum_, const Real comp_)
    : sum(accum_), c(comp_) {}

  Real accum () const { return sum; }
  Real comp () const { return c; }

  KahanBabushkaNeumaierSummer& operator+= (const Real& a) {
    const volatile auto t = sum + a;
    if (std::abs(sum) >= std::abs(a)) {
      const volatile auto z = sum - t;
      const volatile auto y = z + a;
      c += y;
    } else {
      const volatile auto z = a - t;
      const volatile auto y = z + sum;
      c += y;
    }
    sum = t;
    return *this;
  }

  KahanBabushkaNeumaierSummer& operator+= (const KahanBabushkaNeumaierSummer& s) {
    *this += s.sum;
    *this += s.c;
    return *this;
  }

  Real operator() () const { return sum + c; }

  void flip_sign () { sum = -sum; c = -c; }

  Real vimpl () const { return operator()(); }

  static bool unittest(const bool verbose = false);
  
private:
  Real sum, c;
};

struct LongDoubleSummer {
  typedef Real scalar_type;
  typedef long double ld;
  static constexpr int nint = sizeof(ld)/sizeof(int);

  LongDoubleSummer () : sum(0) {}
  LongDoubleSummer (const ld& sum_) : sum(sum_) {}

  LongDoubleSummer& operator+= (const Real& a) {
    sum += static_cast<ld>(a);
    return *this;
  }

  LongDoubleSummer& operator+= (const LongDoubleSummer& s) {
    *this += s.sum;
    return *this;
  }

  Real operator() () const { return static_cast<Real>(sum); }

  void flip_sign () { sum = -sum; }

  const ld& vimpl () const { return sum; }
  
private:
  ld sum;
};

struct OrderInvariantAccumTraits {
  typedef ::Real Real;
  typedef order_invariant_accum::AccumType<Real,true,false>::type Global;
  //todo-pebfb We want <Real,true,true>.
  typedef order_invariant_accum::AccumType<Real,false,true>::type LocalAtomic;
  
  enum : int { global_size_nreal =
               (sizeof(Global) + sizeof(Real) - 1)/sizeof(Real) };
  enum : int { local_atomic_size_nreal =
               (sizeof(LocalAtomic) + sizeof(Real) - 1)/sizeof(Real) };
};

#endif
