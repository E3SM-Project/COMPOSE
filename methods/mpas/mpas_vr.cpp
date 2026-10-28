#include "mpas_vr.hpp"

extern "C" void
remap1c(int nlev, int alg, const double* dp1, const double* dp2, double* Qdpi);

namespace vr {

void homme_vr (int nlev, const Real* dp1, const Real* dp2, Real* Qdpi,
               const bool full_bdy_extrap) {
  const int alg = full_bdy_extrap ? 11 : 10;
  remap1c(nlev, alg, dp1, dp2, Qdpi);
}

static bool assert_grid (const char* name, const Int nb, const Real* b,
                         const bool verbose) {
  bool ok = true;
  for (Int i = 1; i < nb; ++i)
    if (b[i] <= b[i-1]) {
      ok = false;
      if (verbose)
        printf("assert_grid: %s %d %23.16e %23.16e\n",
               name, i, b[i-1], b[i]);
    }
  return ok;
}

bool assert_grids (const Int nb, const Real* bs, const Real* bt,
                   const bool verbose) {
#ifdef NDEBUG
  return true;
#else
  bool ok = true;
  const Real
    tol = (std::abs(bs[nb-1] - bs[0])*
           100*std::numeric_limits<Real>::epsilon()),
    e0 = std::abs(bs[0] - bt[0]),
    e1 = std::abs(bs[nb-1] - bt[nb-1]);
  if (e0 > tol || e1 > tol) {
    ok = false;
    if (verbose)
      printf("assert_grids bounds: %23.16e %23.16e %10.3e, "
             "%23.16e %23.16e %10.3e (%10.3e)\n",
             bs[0], bt[0], e0, bs[nb-1], bt[nb-1], e1, tol);
  }
  if ( ! assert_grid("src", nb, bs, verbose)) ok = false;
  if ( ! assert_grid("tgt", nb, bt, verbose)) ok = false;
  return ok;
#endif
}

} // namespace vr
