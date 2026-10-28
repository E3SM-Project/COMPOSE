#include "mpas_accum.hpp"
#include "mpas_util.hpp"

namespace accum {

static const Int nslot = threaded_bfb_nslot;

static Real accumulate (Real* const a, const Int n, Real slot[nslot]) {
  const Int nperslot = (n + nslot - 1)/nslot;
  ompfor for (Int si = 0; si < nslot; ++si) {
    const Int ib = si*nperslot;
    const Int ie = std::min(ib + nperslot, n);
    Real accum = 0;
    for (Int i = ib; i < ie; ++i) accum += a[i];
    slot[si] = accum;
  }
  reduce_in_place(slot, 1, nslot);
  return slot[0];
}

Real accumulate_threaded_bfb (Real* const a, const Int n, Real slot[nslot],
                              const bool in_threaded_region) {
  Real accum;
  if (in_threaded_region)
    accum = accumulate(a, n, slot);
  else
    omppar accum = accumulate(a, n, slot);
  return accum;
}

// n has to be a power of 2.
void reduce_in_place (Real* a, const Int nscalar, const Int n) {
  Int nit = n >> 1, stride = 1;
  while (nit) {
    const Int stride2 = stride << 1;
    assert(nit*stride2 == n);
#   pragma omp for collapse(2)
    for (Int i = 0; i < n; i += stride2)
      for (int j = 0; j < nscalar; ++j)
        a[i*nscalar + j] += a[(i + stride)*nscalar + j];
    nit >>= 1;
    stride = stride2;
  }
}

} // namespace accum

bool KahanBabushkaNeumaierSummer::unittest (const bool verbose) {
  bool ok = true;
  for (int trial = 0; trial < 6; ++trial) {
    const int sign = trial % 2 == 0 ? 1 : -1;
    const int n = 1 << (trial+10);

    std::vector<Real> v(n);
    for (int i = 0; i < n; ++i) v[i] = sign*std::pow(10*urand(), 1 + 4*urand());
    if (verbose) {
      Real min = 1e30, max = -1e30;
      for (int i = 0; i < n; ++i) max = std::max(max,v[i]);
      for (int i = 0; i < n; ++i) min = std::min(min,v[i]);
      printf("sum: v min,max: %23.16e %23.16e\n", min, max);
    }

    Real sum_naive = 0;
    for (int i = 0; i < n; ++i) sum_naive += v[i];

    for (int test = 0; test < 2; ++test) {
      const int pilim = std::min(50, n);
      std::vector<Real> sums(pilim);
      for (int pi = 0; pi < pilim; ++pi) {
        // Randomly permute v.
        for (int i = 0; i < n; ++i) {
          const int j = urand()*n, k = urand()*n;
          std::swap(v[j], v[k]);
        }

        if (test == 0) {
          // Test += Real.
          KahanBabushkaNeumaierSummer s;
          for (int i = 0; i < n; ++i)
            s += v[i];
          // Test some ctors and getters.
          const auto s1 = s;
          KahanBabushkaNeumaierSummer s2(s1.accum(), s1.comp());
          sums[pi] = s2();
        } else {
          // Also test += KahanBabushkaNeumaierSummer.
          KahanBabushkaNeumaierSummer s1;
          for (int i = 0; i < n/2; ++i)
            s1 += v[i];
          KahanBabushkaNeumaierSummer s2;
          for (int i = n/2; i < n; ++i)
            s2 += v[i];
          s1 += s2;
          sums[pi] = s1();
        }
      }

      static const Real eps = std::numeric_limits<Real>::epsilon();
      if (reldif(sum_naive, sums[0]) > n*eps) {
        ok = false;
        printf("sum FAIL against naive: %6d %23.16e %23.16e %23.16e\n",
               n, sum_naive, sums[0], reldif(sum_naive, sums[0]));
      }
      for (int pi = 1; pi < pilim; ++pi)
        if (reldif(sums[pi], sums[0]) > 2*eps) {
          ok = false;
          printf("sum FAIL against itself: %2d %23.16e %23.16e %23.16e (%23.16e)\n",
                 pi, sum_naive, sums[0], sums[pi], reldif(sums[pi], sums[0]));
        }
      if (verbose)
        printf("sum: %5d %23.16e %23.16e %23.16e\n",
               n, sum_naive, sums[0], reldif(sum_naive, sums[0]));
    }
  }
  return ok;  
}
