#ifndef INCLUDE_MPAS_VR_HPP
#define INCLUDE_MPAS_VR_HPP

#include "mpas_types.hpp"

namespace vr {

// Remap Qdpst from dps to dpt, where dp are the grids' cell widths.
void homme_vr(int nlev, const Real* dps, const Real* dpt, Real* Qdpst,
              const bool full_bdy_extrap = false);

// Assert that the interface values obey the properties needed to use
// homme_vr. Disabled ifdef NDEBUG. Does not call assert but rather returns true
// if OK, false if not.
bool assert_grids(const Int nb, const Real* bs, const Real* bt,
                  const bool verbose=true);

} // namespace vr

#endif
