#ifndef INCLUDE_CEDR_CAAS_INL_HPP
#define INCLUDE_CEDR_CAAS_INL_HPP

namespace cedr {
// ClipAndAssuredSum.
namespace caas {

template <typename ES> KOKKOS_INLINE_FUNCTION
void CAAS<ES>::set_rhom (const Int& lclcellidx, const Real& rhom) {
  d_(lclcellidx) = rhom;
}

template <typename ES> KOKKOS_INLINE_FUNCTION
void CAAS<ES>
::set_Qm (const Int& lclcellidx, const Int& tracer_idx,
          const Real& Qm, const Real& Qm_min, const Real& Qm_max,
          const Real Qm_prev) {
  d_((1 +               tracer_idx)*nlclcells_ + lclcellidx) = Qm;
  d_((1 +   ntracers_ + tracer_idx)*nlclcells_ + lclcellidx) = Qm_min;
  d_((1 + 2*ntracers_ + tracer_idx)*nlclcells_ + lclcellidx) = Qm_max;
}

template <typename ES> KOKKOS_INLINE_FUNCTION
Real CAAS<ES>::get_Qm (const Int& lclcellidx, const Int& tracer_idx) {
  return d_((1 + tracer_idx)*nlclcells_ + lclcellidx);
}

} // namespace caas
} // namespace cedr

#endif
