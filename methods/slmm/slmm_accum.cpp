#include "slmm_defs.hpp"
#include "slmm_accum.hpp"

namespace slmm {

static const Int nslot = accumulate_threaded_bfb_nslot;

static Real accumulate (Real* const a, const Int n, Real slot[nslot]) {
  const Int nperslot = (n + nslot - 1)/nslot;
# pragma omp for
  for (Int si = 0; si < nslot; ++si) {
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
#   pragma omp parallel
    accum = accumulate(a, n, slot);
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

}
