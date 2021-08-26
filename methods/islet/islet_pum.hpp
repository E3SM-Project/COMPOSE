#ifndef INCLUDE_ISLET_PUM_HPP
#define INCLUDE_ISLET_PUM_HPP

#include <map>

#include "islet_types.hpp"
#include "islet_maxeigcomp.hpp"

namespace pum {

struct Options {
  bool threaded;
  Int ne, ntrial, mec_ne;
  Real perturb;
  Options();
};

struct PerturbedUniformMeshMetric : public UserInterpMethod {
  typedef std::shared_ptr<PerturbedUniformMeshMetric> Ptr;

  PerturbedUniformMeshMetric(const InterpMethod::Ptr& im,
                             const Options opts = Options());
  PerturbedUniformMeshMetric(const UserInterpMethod::Ptr& im,
                             const Options opts = Options());
  Real run(Real stop_if_above = 1e3, const bool one_elem_hop_only = false);

  // Can't reset opts.threaded
  void reset_opts(const Options& o);

  // UserInterpMethod interface
  void eval(const Real& x, Real* const v) override;
  const Real* get_xnodes() const override;
  Int get_np() const override;

  // Illustrate why a 1-element hop is the key thing to study.
  void sweep_and_collect_amplitudes(
    const Int npts, const Real threshold,
    // Report meam1 at dx in [0,1] if meam1 >= threshold. This routine does not
    // clear what is already in dx2meam1.
    std::map<Real,Real>& dx2meam1,
    const bool verbose = true);

private:
  Options opts;
  InterpMethod::Ptr base_im;
  MaxEigComputer mec;
  std::vector<std::vector<Real> > xbs, xnodess;

  void init();
};

void demo();

} // namespace pum

#endif
