#ifndef INCLUDE_SLMMIR_D2C_HPP
#define INCLUDE_SLMMIR_D2C_HPP

#include "slmmir.hpp"

// Move between discontinuous and continuous bases for I/O and to support
// mimicking Homme's DSS step. This is thread-scalable, but we don't attempt to
// be efficient: more work is done than necessary in simulating Homme's CG
// basis.
class D2Cer {
  // dglln2cglln_(i) is the CG node associated with DG node i.
  AIdxArray dglln2cglln_;

  // cglln2dglln_(c2dptr_(i) : c2dptr_(i+1)-1) are the DG nodes associated with
  // CG node i.
  AIdxArray cglln2dglln_, c2dptr_;

  // Basis function integrals.
  ARealArray dgbfi_;
  ARealArray cgbfi_;

public:
  typedef std::shared_ptr<D2Cer> Ptr;

  D2Cer(const AIdxArray& dglln2cglln,
        const ARealArray& dgbfi);

  Int get_cnn () const { return len(c2dptr_) - 1; }
  Int get_dnn () const { return len(dglln2cglln_); }
  
  void c2d(const Real* const cg_data, Real* const dg_data) const;
  void d2c(const Real* const dg_data, Real* const cg_data,
           // Inject the value rather than DSS summing. This is useful when
           // working with already continuous data to avoid numerical noise at
           // the eps level from the DSS.
           const bool inject=false) const;
  void dss(Real* const Q, Real* const wrk);

  // Make q continuous despite rho being discontinuous.
  void d2c_q(const Real* const rho_dg, const Real* const q_dg,
             Real* const cg_data);
  void dss_q(const Real* const rho, Real* const q, Real* const wrk);
};

void Qtoq(const Int ntracer, const Int n, const Real* rho, Real* tracers);
void qtoQ(const Int ntracer, const Int n, const Real* rho, Real* tracers);

// Discrete stiffness summation (DSS) to map discontinuous Galerkin field to
// continguous Galerkin.
void dss(D2Cer& d2cer, const bool mixing_ratio_input, const Int ntracer,
         const Int n, Real* wrk, Real* rho, Real* tracers,
         const bool dss_rho = true);

#endif
