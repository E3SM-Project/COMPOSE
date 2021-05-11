#ifndef INCLUDE_SLMMIR_ACCUM_HPP
#define INCLUDE_SLMMIR_ACCUM_HPP

namespace slmm {

static const Int accumulate_threaded_bfb_nslot = 512;

void reduce_in_place(Real* a, const Int nscalar, const Int n);

// fn implements += into its argument.
// slot must be allocated to nslot*nscalar.
template <typename Fn>
void accumulate_threaded_bfb (Fn fn, const Int n, const Int nscalar, Real* slot,
                              Real* accum) {
  static const Int nslot = accumulate_threaded_bfb_nslot;
  const Int nperslot = (n + nslot - 1)/nslot;
# pragma omp parallel for
  for (Int si = 0; si < nslot; ++si) {
    Real* slotp = slot + si*nscalar;
    for (Int j = 0; j < nscalar; ++j) slotp[j] = 0;
    const Int ib = si*nperslot;
    const Int ie = std::min(ib + nperslot, n);
    for (Int i = ib; i < ie; ++i)
      fn(i, slotp);
  }
  reduce_in_place(slot, nscalar, accumulate_threaded_bfb_nslot);
  for (Int j = 0; j < nscalar; ++j) accum[j] = slot[j];
}

template <int nscalar, typename Fn>
void accumulate_threaded_bfb (Fn fn, const Int n, Real* accum) {
  Real slot[accumulate_threaded_bfb_nslot*nscalar];
  accumulate_threaded_bfb(fn, n, nscalar, slot, accum);
}

Real accumulate_threaded_bfb(Real* const a, const Int n,
                             Real slot[accumulate_threaded_bfb_nslot],
                             const bool in_threaded_region = true);

}

#endif
