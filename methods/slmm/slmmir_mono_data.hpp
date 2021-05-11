#ifndef INCLUDE_SLMMIR_MONO_DATA_HPP
#define INCLUDE_SLMMIR_MONO_DATA_HPP

#include <memory>

#include "slmm_gll.hpp"
#include "slmm_gallery.hpp"
using namespace slmm;

#include "slmmir_mesh.hpp"
#include "slmmir_remap_data.hpp"
#include "slmmir.hpp"

// Data for monotonic basis in source cell when such is necessary. Even though
// the cell of interest is a source, and not a target, cell, the calculations
// occur on the Eulerian mesh. Hence MonoData is constant over time iterations.
class MonoData {
  const Filter::Enum filter_type_;
  const Limiter::Enum limiter_type_;
  const Int np_, np2_;
  const AVec3s p_;
  const AIdxs c2n_;
  const ARealArray Ft_;

  Int lqlt_bufsz;
  std::vector<Int> lqlt_tree;
  mutable std::vector<Real> lqlt_buf;

  Int spf_run(const Int n, const Real* w, const Real* a, const Real b,
              const Real* xlo, const Real* xhi, const bool xbds_scalar,
              const Real* y, Real* x, const Limiter::Enum type) const;

public:
  typedef std::shared_ptr<MonoData> Ptr;

  MonoData(const Mesh& m, const RemapData& rd,
           const Filter::Enum monotone_type, const Limiter::Enum limiter_type);

  Filter::Enum global_type () const { return filter_type_; }
  Limiter::Enum local_type () const { return limiter_type_; }

  Int np2 () const { return np2_; }
  const ARealArray dgbfi_mass () const { return Ft_; }

  // Compute q_min, q_max for a source cell.
  void calc_q_min_max(
    const Int& ci, const Real* const rho, const Real* const Q,
    Real& xlo, Real& xhi) const;

  // Positivity limiter.
  Int limit_density(const Int ci,
                    // Data for just cell ti.
                    Real* const rho_tgt, const Real extra_mass=0) const;

  // Compute q_min, q_max for a target cell.
  void calc_q_min_max(
    const RemapData::MT& T, const Int& ti, const Real* const src_rho,
    const Real* const Q_src, Real& q_min, Real& q_max) const;

  // Solve the QP
  //   min_q' sum_i w_i r*_i (q'_i - q*_i)^2
  //    st    sum_i w_i r*_i q'_i = sum_i w_i r*_i q*_i = sum_i w_i Q*_i + extra
  //          q_min <= q'_i <= q_max.
  Int limit_tracer(const Int ti, Real q_min, Real q_max,
                   // Data for just cell ti.
                   const Real* const tgt_rho, Real* const Q_tgt,
                   const bool expand_bounds=false, const Real Q_extra=0,
                   const bool nochecks=false,
                   Limiter::Enum limiter_type = Limiter::none) const;

  Int limit_tracer(
    const RemapData::MT& T, const Int ti,
    const Real* const src_rho, const Real* const Q_src,
    const Real* const tgt_rho, Real* const Q_tgt,
    const bool expand_bounds=false, const Real Q_extra=0,
    const bool nochecks=false) const;

  Int limit_tracer(const Int ti, Real* q_min, Real* q_max,
                   // Data for just cell ti.
                   const Real* const tgt_rho, Real* const Q_tgt,
                   const bool expand_bounds_allowed, const Real Qm_extra,
                   const bool nochecks=false) const;
};

#endif
