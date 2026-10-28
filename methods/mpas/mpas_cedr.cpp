#include "mpas_cedr.hpp"
#include "mpas_comm_util.hpp"

void caas (Comm& he, const Int nprob, const Int m, const Int n,
           const Real* const a, const bool is_3d, const Acr2& dzs,
           const Acr2& zt, const Acr3& xlo, const Acr3& xhi,
           const Acr3& xmassref, const Acr3& y, const Ar3& x, const Ar1& wrk,
           const bool clip, const bool bfb) {
  static const int nsum = 3, bos = nsum-2;
  const auto f = [&] (const int ci, Real* acc) {
    for (Int i = 0; i < nprob*nsum; ++i)
      acc[i] = 0;
    for (int qi = 0; qi < nprob; ++qi) {
      const int os = qi*nsum;
      for (Int k = 0; k < m; ++k) {
        const Real dzski = is_3d ? dzs(k,ci) : 1;
        const Real dztki = is_3d ? zt(k+1,ci) - zt(k,ci) : 1;
        x(qi,k,ci) = std::max(xlo(qi,k,ci), std::min(xhi(qi,k,ci), y(qi,k,ci)));
        const Real xi = x(qi,k,ci), ai = a[ci];
        acc[os + 0    ] += ai*(dzski*xmassref(qi,k,ci) - dztki*xi);
        acc[os + bos  ] += ai*dztki*(xhi(qi,k,ci) - xi);
        acc[os + bos+1] += ai*dztki*(xi - xlo(qi,k,ci));
      }
    }
  };
  Real* sum = wrk.data() + nsum*nprob;
  repro_sum(he, n, nsum*nprob, wrk, f, bfb, sum);
  for (Int qi = 0; qi < nprob; ++qi) {
    const Int os = nsum*qi;
    const Real dm = sum[os];
    if (dm > 0) {
      if (sum[os+bos] == 0) continue;
      const Real f = dm/sum[os+bos];
      ompparfor for (Int ci = 0; ci < n; ++ci)
        for (Int k = 0; k < m; ++k)
          x(qi,k,ci) += f*(xhi(qi,k,ci) - x(qi,k,ci));
    } else {
      if (sum[os+bos+1] == 0) continue;
      const Real f = dm/sum[os+bos+1];
      ompparfor for (Int ci = 0; ci < n; ++ci)
        for (Int k = 0; k < m; ++k)
          x(qi,k,ci) += f*(x(qi,k,ci) - xlo(qi,k,ci));
    }
  }
  // Clip again for numerics.
  if (clip)
    ompparfor for (Int ci = 0; ci < n; ++ci)
      for (Int k = 0; k < m; ++k)
        for (Int qi = 0; qi < nprob; ++qi)
          x(qi,k,ci) = std::max(xlo(qi,k,ci), std::min(xhi(qi,k,ci), x(qi,k,ci)));
}
