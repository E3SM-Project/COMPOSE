#include "slmmir_mono_data.hpp"

MonoData::MonoData (const Mesh& m, const RemapData& rd,
                    const Filter::Enum monotone_type,
                    const Limiter::Enum limiter_type)
  : filter_type_(monotone_type), limiter_type_(limiter_type), np_(m.np),
    np2_(square(np_)), p_(m.geo_p), c2n_(m.geo_c2n), Ft_(rd.dgbfi_mass())
{
  const auto nelem = nslices(m.geo_c2n);
  Array<int> rowptr(nelem + 1), colidx(nelem);
  for (int i = 0; i < nelem; ++i) {
    rowptr[i] = i;
    colidx[i] = i;
  }
  rowptr[nelem] = nelem;
  if (limiter_type == Limiter::qlt) {
    lqlt_bufsz = spf::local_qlt_tensor2d_init(np_, lqlt_tree);
    lqlt_buf.resize(omp_get_max_threads() * lqlt_bufsz);
  }
}

// Compute q_min, q_max for a source cell.
void MonoData
::calc_q_min_max (
  const Int& ci, const Real* const rho, const Real* const Q,
  Real& xlo, Real& xhi) const
{
  xlo = 1;
  xhi = 0;
  for (Int i = 0; i < np2_; ++i) {
    const Real q = rho[i] <= 0 ? 1 : Q[i]/rho[i];
    xlo = std::min(xlo, q);
    xhi = std::max(xhi, q);
  }
#ifdef EXPENSIVE_CHECKS
  //todo Figure out if these hold in FP or if there is a 1 or 2 ulp error
  // possible.
  assert(xlo >= 0);
  if (xhi > 1) prc(xhi-1);
  assert(xhi <= 1);
#endif
}

Int MonoData
::spf_run (const Int n, const Real* w, const Real* a, const Real b,
           const Real* xlo, const Real* xhi, const bool xbds_scalar,
           const Real* y, Real* x, const Limiter::Enum type) const {
  switch (type) {
  case Limiter::mn2: return spf::solve_1eq_bc_qp(n, w, a, b, xlo, xhi, xbds_scalar, y, x);
  case Limiter::caas: return spf::clip_and_sum(n, w, a, b, xlo, xhi, xbds_scalar, y, x);
  case Limiter::caags: return spf::clip_and_weighted_sum(n, w, a, b, xlo, xhi, xbds_scalar, y, x);
  case Limiter::qlt: {
    SIQK_THROW_IF(lqlt_tree.empty() || lqlt_buf.empty(),
                  "Local QLT requires initialization.");
    return spf::local_qlt_tensor2d_run(lqlt_tree.data(),
                                       lqlt_buf.data() + omp_get_thread_num() * lqlt_bufsz,
                                       w, a, b, xlo, xhi, xbds_scalar, y, x);
  }
  default:
    assert(0);
    return -1;
  }
}

// Positivity limiter.
Int MonoData
::limit_density (const Int ci,
                 // Data for just cell ti.
                 Real* const rho_tgt, const Real extra_mass) const {
  const Int ns = np2_*ci;
  Real mass_tgt = extra_mass;
  bool any_below = false;
  for (Int i = 0; i < np2_; ++i) {
    if (rho_tgt[i] < 0) any_below = true;
    mass_tgt += rho_tgt[i] * Ft_[ns+i];
  }
  Int info = mass_tgt < 0 ? -1 : 0;
  if ( ! any_below && extra_mass == 0) return info;
  Real mass = 0;
  for (Int i = 0; i < np2_; ++i) {
    if (rho_tgt[i] < 0) rho_tgt[i] = 0;
    mass += rho_tgt[i] * Ft_[ns+i];
  }
  Real delta_mass = mass_tgt - mass;
  if (delta_mass >= 0) {
    // This is linearly invariant.
    Real fac = 0;
    for (Int i = 0; i < np2_; ++i)
      fac += Ft_[ns+i];
    fac = delta_mass/fac;
    for (Int i = 0; i < np2_; ++i)
      rho_tgt[i] += fac;
  } else {
    // This case cannot be, since 0 is absolute. Adding a constant b
    // everywhere will give results variant to b. However, by using
    // Limiter::mn2, it is invariant if no clipping occurs. (CAAGS with
    // appropriate user weight would also be.)
    Real zeros[GLL::np_max*GLL::np_max] = {0};
    Real y[GLL::np_max*GLL::np_max];
    for (Int i = 0; i < np2_; ++i) y[i] = rho_tgt[i];
    info = spf_run(np2_, &Ft_[ns], &Ft_[ns], mass_tgt, zeros, y, false, y,
                   rho_tgt, Limiter::mn2);
  }
  return info;
}

// Compute q_min, q_max for a target cell.
void MonoData::calc_q_min_max (
  const RemapData::MT& T, const Int& ti, const Real* const src_rho,
  const Real* const Q_src, Real& q_min, Real& q_max) const
{
  const Size* const rowptr = T.rowptr();
  const Int* const colidx = T.colidx();
  q_min = 1;
  q_max = 0;
  for (Size j = rowptr[ti]; j < rowptr[ti+1]; ++j) {
    const Int ti_src = colidx[j];
    for (Int i = 0, k = np2_*ti_src; i < np2_; ++i, ++k) {
      const Real rho = src_rho[k];
      const Real q_src = rho <= 0 ? 1 : Q_src[k] / rho;
      q_min = std::min(q_min, q_src);
      q_max = std::max(q_max, q_src);
    }
  }
  // q_min should be >= 0 and q_max <= 1.
  q_min = std::max(q_min, 0.0);
  q_max = std::min(q_max, 1.0);
}

// Solve the QP
//   min_q' sum_i w_i r*_i (q'_i - q*_i)^2
//    st    sum_i w_i r*_i q'_i = sum_i w_i r*_i q*_i = sum_i w_i Q*_i + extra
//          q_min <= q'_i <= q_max.
Int MonoData
::limit_tracer (const Int ti, Real q_min, Real q_max,
                // Data for just cell ti.
                const Real* const tgt_rho, Real* const Q_tgt,
                const bool expand_bounds, const Real Q_extra,
                const bool nochecks, Limiter::Enum limiter_type) const
{
  const Int ns = np2_*ti;
  if (limiter_type == Limiter::none) limiter_type = limiter_type_;

  // Compute tracer and total masses in the target cell. (The total mass is
  // computed redundantly for each tracer; might change that.)
  Real tracer_mass_tgt = 0, total_mass_tgt = 0;
  for (Int i = 0; i < np2_; ++i) {
    total_mass_tgt += tgt_rho[i] * Ft_[ns+i];
    tracer_mass_tgt += Q_tgt[i] * Ft_[ns+i];
  }
  tracer_mass_tgt += Q_extra;

  if ( ! nochecks || expand_bounds) {
    const Real Qmorm = tracer_mass_tgt / total_mass_tgt;
    if (expand_bounds) {
      q_min = std::min(q_min, Qmorm);
      q_max = std::max(q_max, Qmorm);
    }
    // Cell mean boundedness is necessary and sufficient for QP feasibility.
    if ( ! expand_bounds && (tracer_mass_tgt < q_min*total_mass_tgt ||
                             tracer_mass_tgt > q_max*total_mass_tgt)) {
      // A constant field is a problem because the cell mean will almost surely
      // be out of bounds in the last digit. Treat it specially.
      const bool possible_constant_field = q_min == q_max;
      if ( ! ( // The target field could be constant because the bounds are the
             // same, indicating the source field is constant.
             possible_constant_field &&
             // If <Q>/<rho> is sufficiently close, consider it OK.
             (std::abs(Qmorm - q_min) <=
              10*std::numeric_limits<Real>::epsilon()))) {
        // The necessary and sufficient condition really is not met, so return
        // the failure status.
        //pr(puf(tracer_mass_tgt) pu(Q_extra) pu(q_min) pu(q_max-q_min) pu(Qmorm));
        return -2;
      }
    }
  }

  Int info;
  { // Set up and solve the QP.
    Real a[GLL::np_max*GLL::np_max], y[GLL::np_max*GLL::np_max],
      x[GLL::np_max*GLL::np_max];
    for (Int i = 0; i < np2_; ++i) a[i] = Ft_[ns+i] * tgt_rho[i];
#ifndef QP_AW
    const Real* const w = a;
#else
    Real w[GLL::np_max*GLL::np_max];
    for (Int i = 0; i < np2_; ++i) w[i] = Ft_[ns+i] * square(tgt_rho[i]);
#endif
    for (Int i = 0; i < np2_; ++i) y[i] = Q_tgt[i] / tgt_rho[i];
    const Real b = tracer_mass_tgt;
    const Real xlo = q_min, xhi = q_max;
    copy(x, y, np2_);
    info = spf_run(np2_, w, a, b, &xlo, &xhi, true, y, x, limiter_type);
    if (expand_bounds) info = 1;
#ifdef EXPENSIVE_CHECKS
    {
      const bool ok =
        spf::check_1eq_bc_qp_foc("limit_tracer", np2_, w, a, b, &xlo, &xhi,
                                 true, y, x, false);
      if (info < 0 || ! ok) {
        prc(info); prc(xlo); prc(xhi);
        prarr("w", w, np2_); prarr("y", y, np2_); prarr("x", x, np2_);
        prarr("tgt_rho", tgt_rho, np2_);
        Real num = 0, den = 0, blo = 0, bhi = 0;
        for (Int i = 0; i < np2_; ++i) {
          num += square(y[i] - x[i]);
          den += square(y[i]);
          blo += a[i] * xlo;
          bhi += a[i] * xhi;
        }
        prc(b); prc(blo); prc(bhi); prc(b-blo); prc(bhi-b);
        pr("reldif(y,x) | " << std::sqrt(num/den));
      }
    }
#endif
#if 0
    // For safety against ~1 ulp errors, since 0 and 1 are special, put q in
    // [0,1] bounds.
    for (Int i = 0; i < np2_; ++i) {
      x[i] = std::min<Real>(1, x[i]);
      x[i] = std::max<Real>(0, x[i]);
    }
#endif
    if (nochecks || info == 1)
      for (Int i = 0; i < np2_; ++i)
        Q_tgt[i] = x[i]*tgt_rho[i];
  }
  return info;
}

Int MonoData::limit_tracer (
  const RemapData::MT& T, const Int ti,
  const Real* const src_rho, const Real* const Q_src,
  const Real* const tgt_rho, Real* const Q_tgt,
  const bool expand_bounds, const Real Q_extra, const bool nochecks) const
{
  Real q_min, q_max;
  calc_q_min_max(T, ti, src_rho, Q_src, q_min, q_max);
  return limit_tracer(ti, q_min, q_max, tgt_rho + ti*np2_, Q_tgt + ti*np2_,
                      expand_bounds, Q_extra, nochecks);
}

static void
expand_bounds (const Int N, Real* q_min, Real* q_max,
               const Real* wgt, const Real* rho, Real Qm_extra)
{
  Real* const q_bnd = Qm_extra < 0 ? q_min : q_max;
  // Try solving a QP that adjusts a q bound.
  Real Qm = Qm_extra;
  Real w[GLL::np_max*GLL::np_max], q_bnd_min[GLL::np_max*GLL::np_max],
    q_bnd_max[GLL::np_max*GLL::np_max], q_bnd_orig[GLL::np_max*GLL::np_max];
  q_bnd_min[0] = q_min[0];
  q_bnd_max[0] = q_max[0];
  for (Int i = 0; i < N; ++i) {
    const Real rhom = wgt[i]*rho[i];
    Qm += q_bnd[i]*rhom;
    q_bnd_orig[i] = q_bnd[i];
    w[i] = rhom;
    if (Qm_extra < 0) {
      q_bnd_min[0] = std::min(q_bnd_min[0], q_min[i]);
      q_bnd_max[i] = q_max[i];
    } else {
      q_bnd_min[i] = q_min[i];
      q_bnd_max[0] = std::max(q_bnd_max[0], q_max[i]);
    }
  }
  if (Qm_extra < 0)
    for (Int i = 1; i < N; ++i) q_bnd_min[i] = q_bnd_min[0];
  else
    for (Int i = 1; i < N; ++i) q_bnd_max[i] = q_bnd_max[0];
  // Check for feasibility.
  bool feasible; {
    Real Qm_lo = 0, Qm_hi = 0;
    for (Int i = 0; i < N; ++i) {
      Qm_lo += q_bnd_min[i]*w[i];
      Qm_hi += q_bnd_max[i]*w[i];
    }
    feasible = Qm_lo <= Qm && Qm <= Qm_hi;
  }
  if (feasible) {
    slmm::spf::solve_1eq_bc_qp(N, w, w, Qm, q_bnd_min, q_bnd_max, false,
                               q_bnd_orig, q_bnd);
  } else {
    // The QP isn't feasible, so set the bound to the uniform value required
    // to get feasibility.
    Real rhom_tot = 0, Qm_tot = Qm_extra;
    for (Int i = 0; i < N; ++i) {
      const Real rhom = wgt[i]*rho[i];
      rhom_tot += rhom;
      Qm_tot += q_bnd_orig[i]*rhom;
    }
    const Real q_tot = Qm_tot / rhom_tot;
    for (Int i = 0; i < N; ++i)
      q_bnd[i] = q_tot;
  }
}

Int MonoData
::limit_tracer (const Int ti, Real* q_min, Real* q_max,
                // Data for just cell ti.
                const Real* const tgt_rho, Real* const Q_tgt,
                const bool expand_bounds_allowed, const Real Qm_extra,
                const bool nochecks) const
{
  const Int ns = np2_*ti;

  Real Qm_tot = Qm_extra, rhom_tot = 0, Qm_min = 0, Qm_max = 0,
    q_min_all = 1, q_max_all = 0;
  for (Int i = 0; i < np2_; ++i) {
    const Real rhom = tgt_rho[i] * Ft_[ns+i];
    rhom_tot += rhom;
    Qm_tot += Q_tgt[i] * Ft_[ns+i];
    Qm_min += q_min[i] * rhom;
    Qm_max += q_max[i] * rhom;
    q_min_all = std::min(q_min_all, q_min[i]);
    q_max_all = std::max(q_max_all, q_max[i]);
  }

  if ( ! nochecks || expand_bounds_allowed) {
    const bool lo = Qm_tot < Qm_min, hi = Qm_tot > Qm_max;
    if (lo || hi) {
      if ( ! expand_bounds_allowed) return -1;
      expand_bounds(np2_, q_min, q_max, &Ft_[ns], tgt_rho,
                    Qm_tot - (lo ? Qm_min : Qm_max));
    }
  }

  Int info;
  { // Set up and solve the QP.
    Real a[GLL::np_max*GLL::np_max], y[GLL::np_max*GLL::np_max],
      x[GLL::np_max*GLL::np_max];
    for (Int i = 0; i < np2_; ++i) a[i] = Ft_[ns+i] * tgt_rho[i];
#ifndef QP_AW
    const Real* const w = a;
#else
    Real w[GLL::np_max*GLL::np_max];
    for (Int i = 0; i < np2_; ++i) w[i] = Ft_[ns+i] * square(tgt_rho[i]);
#endif
    for (Int i = 0; i < np2_; ++i) y[i] = Q_tgt[i] / tgt_rho[i];
    const Real b = Qm_tot;
    copy(x, y, np2_);
    info = spf_run(np2_, w, a, b, q_min, q_max, false, y, x, limiter_type_);
    if (info < 0 && expand_bounds_allowed) {
      // FP must be the problem.
      info = 1;
    }
    // For safety against ~1 ulp errors, since 0 and 1 are special, put q in
    // [0,1] bounds.
#if 0
    for (Int i = 0; i < np2_; ++i) {
      x[i] = std::min<Real>(1, x[i]);
      x[i] = std::max<Real>(0, x[i]);
    }
#endif
    if (nochecks || info == 1)
      for (Int i = 0; i < np2_; ++i)
        Q_tgt[i] = x[i]*tgt_rho[i];
  }
  return info;
}
