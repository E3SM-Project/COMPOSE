#ifndef INCLUDE_MPAS_CEDR_HPP
#define INCLUDE_MPAS_CEDR_HPP

#include "mpas_types.hpp"
#include "mpas_comm.hpp"

/* Find a solution to
     dzt := zt(1:m) - zt(0:m-1)
     min_x norm(dzt x(qi,:,:) - dzt y(qi,:,:), 1)
      st   sum_{k=0:m-1, i=0:n-1}
             a[i] dzt(k,i) x(qi,k,i)
               = sum_{i,k} a[i] dzs(k,i) xmassref(qi,k,i)
           xlo(qi,:,:) <= x(qi,:,:) <= xhi(qi,:,:)
   for qi in 0:nprob-1.
 */
void caas(Comm& he, const Int nprob, const Int m, const Int n,
          const Real* const a, // a[0:n-1]
          // dzs is (0:m-1, 0:n-1); dzt is (0:m, 0:n-1).
          const bool is_3d, const Acr2& dzs, const Acr2& zt,
          // 3D arrays have dims (0:nprob-1, 0:m-1, 0:n-1).
          const Acr3& xlo, const Acr3& xhi,
          const Acr3& xmassref, const Acr3& y, const Ar3& x,
          // If !bfb, wrk has length at least 6*nprob.
          // If  bfb, wrk has length at least 3*nprob*(n+nthr+1), where nthr is
          //          the number of threads to use.
          const Ar1& wrk, const bool clip, const bool bfb);

#endif
