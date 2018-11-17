#ifndef INCLUDE_CEDR_TEST_RANDOMIZED_INL_HPP
#define INCLUDE_CEDR_TEST_RANDOMIZED_INL_HPP

#include "cedr_test_randomized.hpp"

namespace cedr {
namespace test {

template <typename CDRT, typename ExeSpace>
Int TestRandomized::run (const Int nrepeat, const bool write) {
  const Int nt = tracers_.size(), nlclcells = gcis_.size();

  Values v(nt, nlclcells);
  generate_rho(v);
  for (const auto& t : tracers_) {
    generate_Q(t, v);
    perturb_Q(t, v);
  }

  if (write)
    for (const auto& t : tracers_)
      write_pre(t, v);

  CDRT cdr = static_cast<CDRT&>(get_cdr());
  ValuesDevice<ExeSpace> vd(v);
  const auto set_rhom = KOKKOS_LAMBDA (const Int& i) {
    cdr.set_rhom(i, 0, vd.rhom()[i]);
  };
  Kokkos::parallel_for(nlclcells, set_rhom);
  for (Int trial = 0; trial <= nrepeat; ++trial) {
    for (Int ti = 0; ti < nt; ++ti) {
      Real* Qm_min = v.Qm_min(ti), * Qm = v.Qm(ti), * Qm_max = v.Qm_max(ti),
        * Qm_prev = v.Qm_prev(ti);
      for (Int i = 0; i < nlclcells; ++i)
        cdr.set_Qm(i, ti, Qm[i], Qm_min[i], Qm_max[i], Qm_prev[i]);
    }

    run_impl(trial);
  }

  for (Int ti = 0; ti < nt; ++ti) {
    Real* Qm = v.Qm(ti);
    for (Int i = 0; i < nlclcells; ++i)
      Qm[i] = cdr.get_Qm(i, ti);
  }

  if (write)
    for (const auto& t : tracers_)
      write_post(t, v);
  return check(cdr_name_, *p_, tracers_, v);
}

} // namespace test
} // namespace cedr

#endif
