#include "slmmir_snapshot.hpp"
#include "slmmir_p_refine.hpp"
#include "slmm_array.hpp"
#include "slmm_accum.hpp"

void check (const Snapshot& s, const Mesh& m, const Real* q, const Real* w, const Int dnn) {
  const Int ne = dnn/square(m.np), nic = s.q.size(), nscalar = 4*nic;
  const Int np = m.np, np2 = square(np), np_ref = s.m->np, np2_ref = square(np_ref);
  Interpolator2D::Ptr i2d;
  if (np_ref != np)
    i2d = std::make_shared<Interpolator2D>(np_ref, s.m->basis, m.np, m.basis);
  Array<Real> slot((accumulate_threaded_bfb_nslot+1)*nscalar);
  const auto err_fn = [&] (const Int ie, Real* accum) {
    for (Int iq = 0; iq < nic; ++iq) {
      const Real* qi_ref = &s.q[iq][np2_ref*ie];
      const Real* qi_ref_i2d = qi_ref;
      Real buf[Basis::np_max*Basis::np_max];
      if (i2d) {
        qi_ref_i2d = buf;
        i2d->apply(qi_ref, buf); 
      }
      const Real* qi = q + iq*dnn + ie*np2;
      Real* const a = accum + 4*iq;
      for (Int i = 0; i < np2; ++i) {
        const Int k = np2*ie + i;
        const auto d = qi[i] - qi_ref_i2d[i];
        a[0] += w[k]*std::abs(d);
        a[1] += w[k]*std::abs(qi_ref_i2d[i]);
        a[2] += w[k]*square(d);
        a[3] += w[k]*square(qi_ref_i2d[i]);
      }
    }
  };
  Real* err = slot.data();
  accumulate_threaded_bfb(err_fn, ne, nscalar, slot.data() + nscalar, err);
  printf("\n");
  for (Int iq = 0; iq < nic; ++iq) {
    Real* const e = err + 4*iq;
    printf("M %.3s |                              %1.3e %1.3e\n",
           gallery::InitialCondition::to_string(s.ics[iq]),
           e[0]/e[1], std::sqrt(e[2]/e[3]));
  }
}
