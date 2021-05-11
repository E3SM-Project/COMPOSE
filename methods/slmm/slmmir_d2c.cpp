#include "slmmir_d2c.hpp"

void Qtoq (const Int nics, const Int n, const Real* rho, Real* tracers) {
  for (Int iic = 0; iic < nics; ++iic) {
    Real* const t = tracers + iic*n;
#   pragma omp parallel for
    for (Int i = 0; i < n; ++i)
      t[i] /= rho[i];
  }
}

void qtoQ (const Int nics, const Int n, const Real* rho, Real* tracers) {
  for (Int iic = 0; iic < nics; ++iic) {
    Real* const t = tracers + iic*n;
#   pragma omp parallel for
    for (Int i = 0; i < n; ++i)
      t[i] *= rho[i];
  }
}

D2Cer::D2Cer (const AIdxArray& dglln2cglln,
              const ARealArray& dgbfi)
  : dglln2cglln_(dglln2cglln), dgbfi_(dgbfi)
{
  Int cnn = 0;
  for (Int i = 0; i < len(dglln2cglln); ++i)
    cnn = std::max(cnn, dglln2cglln(i));
  ++cnn;

  resize(cgbfi_, cnn);
  copy(cgbfi_, 0);
  for (Int i = 0; i < len(dglln2cglln); ++i)
    cgbfi_(dglln2cglln(i)) += dgbfi_(i);

  // Like setting up a CSR matrix transpose.
  resize(cglln2dglln_, len(dglln2cglln));
  resize(c2dptr_, cnn+1);
  for (Int i = 0; i <= cnn; ++i)
    c2dptr_(i) = 0;
  // Count.
  for (Int i = 0; i < len(dglln2cglln); ++i) {
    const Int i_cgll = dglln2cglln(i);
    if (i_cgll < cnn-1)
      ++c2dptr_(i_cgll+2);
  }
  // Cumsum.
  for (Int i = 3; i <= cnn; ++i)
    c2dptr_(i) += c2dptr_(i-1);
  // Populate.
  for (Int i = 0; i < len(dglln2cglln); ++i)
    cglln2dglln_(c2dptr_(dglln2cglln(i)+1)++) = i;
}

void D2Cer::c2d (const Real* const cg_data, Real* const dg_data) const {
# pragma omp parallel for
  for (Int i = 0; i < len(dglln2cglln_); ++i)
    dg_data[i] = cg_data[dglln2cglln_(i)];
}

void D2Cer::d2c (const Real* const dg_data, Real* const cg_data,
                 const bool inject) const {
  const Int cnn = len(c2dptr_) - 1;
# pragma omp parallel for
  for (Int ic = 0; ic < cnn; ++ic) {
    if (inject) {
      cg_data[ic] = dg_data[cglln2dglln_(c2dptr_(ic))];
    } else {
      Real num = 0, den = 0;
      const Int id = cglln2dglln_(ic);
      Real rho_min = dg_data[id], rho_max = rho_min;
      for (Int j = c2dptr_(ic); j < c2dptr_(ic+1); ++j) {
        const Int id = cglln2dglln_(j);
        const auto rho = dg_data[id];
        num += dgbfi_(id) * rho;
        den += dgbfi_(id);
        rho_min = std::min(rho_min, rho);
        rho_max = std::max(rho_max, rho);
      }
      cg_data[ic] = num / den;
      // See comment below regarding this step.
      cg_data[ic] = std::min(rho_max, std::max(rho_min, cg_data[ic]));
    }
  }
}

void D2Cer::dss (Real* const Q, Real* const wrk) {
  d2c(Q, wrk);
  c2d(wrk, Q);
}

void D2Cer::d2c_q (const Real* const rho_dg, const Real* const q_dg,
                   Real* const cg_data) {
  const Int cnn = len(c2dptr_) - 1;
# pragma omp parallel for
  for (Int ic = 0; ic < cnn; ++ic) {
    Real num = 0, den = 0;
    const Int id = cglln2dglln_(ic);
    Real q_min = q_dg[id], q_max = q_min;
    for (Int j = c2dptr_(ic); j < c2dptr_(ic+1); ++j) {
      const Int id = cglln2dglln_(j);
      const auto q = q_dg[id];
      num += dgbfi_(id) * rho_dg[id] * q;
      den += dgbfi_(id) * rho_dg[id];
      q_min = std::min(q_min, q);
      q_max = std::max(q_max, q);
    }
    cg_data[ic] = num / den;
    // A weighted sum can't give new extremal values, so any new such value is
    // caused by numerical error. Limit this. The trade off is that mass
    // conservation takes a slight hit, at a level that can't be fixed by
    // sequencing the global limiter to integrate source mass before the DSS. In
    // HOMME, I would recommend not clipping the DSS for numerical noise.
    cg_data[ic] = std::min(q_max, std::max(q_min, cg_data[ic]));
  }
}

void D2Cer::dss_q (const Real* const rho, Real* const q, Real* const wrk) {
  d2c_q(rho, q, wrk);
  c2d(wrk, q);  
}

void dss (D2Cer& d2cer, const bool mixing_ratio_input, const Int ntracer,
          const Int n, Real* wrk, Real* rho, Real* tracers, const bool dss_rho) {
  if (mixing_ratio_input) {
    if (dss_rho) {
      qtoQ(ntracer, n, rho, tracers);
      d2cer.dss(rho, wrk);
      for (Int ti = 0; ti < ntracer; ++ti)
        d2cer.dss(tracers + ti*n, wrk);
      Qtoq(ntracer, n, rho, tracers);
    } else {
      for (Int ti = 0; ti < ntracer; ++ti)
        d2cer.dss_q(rho, tracers + ti*n, wrk);
    }
  } else {
    if (dss_rho) d2cer.dss(rho, wrk);
    for (Int ti = 0; ti < ntracer; ++ti)
      d2cer.dss(tracers + ti*n, wrk);
  }
}
