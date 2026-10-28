#include "mpas_comm.hpp"

void Comm
::calc_maxmag (const int nterm, const int nscalar, const Real* const v,
               Real* maxmag_) {
  static constexpr int max_nscalar = 128; // for OpenMP reduction
  assert(nscalar <= max_nscalar);
  Real maxmag[max_nscalar] = {0};
#if defined _OPENMP && _OPENMP >= 201500
  // Need OpenMP 4.5 for array.
# pragma omp parallel reduction (max: maxmag)
#endif
  for (int i = 0; i < nterm; ++i)
    for (int j = 0; j < nscalar; ++j)
      maxmag[j] = std::max(maxmag[j], std::abs(v[i*nscalar + j]));
  for (int j = 0; j < nscalar; ++j)
    maxmag_[j] = maxmag[j];
}
