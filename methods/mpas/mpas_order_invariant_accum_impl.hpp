#ifndef INCLUDE_MPAS_ORDER_INVARIANT_ACCUM_IMPL_HPP
#define INCLUDE_MPAS_ORDER_INVARIANT_ACCUM_IMPL_HPP

#include <type_traits>

#include "mpas_order_invariant_accum.hpp"

namespace order_invariant_accum {

template <typename IRepT>
Accumulator<IRepT> Accumulator<IRepT>
::sum (const Ar1& wrk, const int nterm, const Real* terms, Real maxmag) {
  // If maxmag isn't provided, determine it.
  if (maxmag == -1) {
#ifdef _OPENMP
#   pragma omp parallel reduction (max: maxmag)
#endif
    for (int i = 0; i < nterm; ++i)
      maxmag = std::max(maxmag, std::abs(terms[i]));
    maxmag *= nterm;
  }
  // Dispatch a team of accumulators, each of which sums a portion of the
  // terms serially.
  const int nthr = mpas_omp_get_max_threads();
  const int naccs = std::min(nthr, int((wrk.size()*sizeof(Real))/sizeof(IRep)));
  assert(naccs > 0);
  auto* const ireps = reinterpret_cast<IRep*>(wrk.data());
  omppar {
    const int tid = mpas_omp_get_thread_num();
    if (tid < naccs) {
      const int n = (nterm + naccs - 1)/naccs;
      const int ibeg = tid*n, iend = std::min(nterm, (tid+1)*n);
      Accumulator slcl(maxmag);
      for (int i = ibeg; i < iend; ++i)
        slcl += terms[i];
      ireps[tid] = slcl.get_irep();
    }
  }
  // Sum accumulators in serial.
  Accumulator stot(maxmag);
  for (int i = 0; i < naccs; ++i)
    stot += Accumulator(maxmag, ireps[i]);
    
  return stot;
}

template <typename IRepT>
void Accumulator<IRepT>
::sum (const Ar1& wrk, const Real* const terms, const Int nterm,
       Accumulator<IRepT>* accum, const Int nscalar,
       const bool assert_on_nthr_too_small) {
  const int nthr = mpas_omp_get_max_threads();
  const int nuse = std::min(nthr, int((wrk.size()*sizeof(Real))/sizeof(IRep)));
  assert(nuse > 0);
  if (assert_on_nthr_too_small) assert(nuse == nthr);
  auto* const slcls = reinterpret_cast<Accumulator*>(wrk.data());
  omppar {
    const int tid = mpas_omp_get_thread_num();
    if (tid < nuse) {
      const int n = (nterm + nuse - 1)/nuse;
      const int ibeg = tid*n, iend = std::min(nterm, (tid+1)*n);
      auto* const slcl = slcls + tid*nscalar;
      for (int j = 0; j < nscalar; ++j) {
        slcl[j].maxmag = accum[j].get_maxmag();
        slcl[j].irep = 0;
      }
      for (int i = ibeg; i < iend; ++i) {
        for (int j = 0; j < nscalar; ++j)
          slcl[j] += terms[i*nscalar + j];
      }
    }
  }
  // Sum accumulators in serial.
  for (int i = 0; i < nuse; ++i)
    for (int j = 0; j < nscalar; ++j)
      accum[j] += slcls[i*nscalar + j];
}

} // namespace order_invariant_accum

#endif
