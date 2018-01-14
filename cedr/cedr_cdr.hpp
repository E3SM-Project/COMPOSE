#ifndef INCLUDE_CEDR_CDR_HPP
#define INCLUDE_CEDR_CDR_HPP

#include "cedr_mpi.hpp"

namespace cedr {
// Constrained Density Reconstructor interface.
struct CDR {
  // Set up QLT tracer metadata. Once end_tracer_declarations is called, it is
  // an error to call declare_tracer again. Call declare_tracer in order of the
  // tracer index in the caller's numbering. It is an error to call this
  // function from a parallel region.
  virtual void declare_tracer(int problem_type) = 0;

  // It is an error to call this function from a parallel region.
  virtual void end_tracer_declarations() = 0;

  virtual int get_problem_type(const Int& tracer_idx) const = 0;

  virtual Int get_num_tracers() const = 0;

  // set_{rhom,Qm}: Set cell values prior to running the QLT algorithm.
  //   set_rhom must be called before set_Qm.
  //   Notation:
  //     rho: Total density.
  //       Q: Tracer density.
  //       q: Tracer mixing ratio = Q/rho.
  //      *m: Mass corresponding to the density; results from an integral over a
  //          region, such as a cell.
  //   Some CDRs have a nontrivial local <-> global cell index map. For these
  // CDRs, lclcellidx may be nontrivial. For others, the caller should provide
  // the index into the local cell.
  virtual void set_rhom(
    const Int& lclcellidx,
    // Current total mass in this cell.
    const Real& rhom) = 0;

  virtual void set_Qm(
    const Int& lclcellidx, const Int& tracer_idx,
    // Current tracer mass in this cell.
    const Real& Qm,
    // Minimum and maximum permitted tracer mass in this cell.
    const Real& Qm_min, const Real& Qm_max,
    // If mass conservation is requested, provide the previous Qm, which will be
    // summed to give the desired global mass.
    const Real Qm_prev = -1) = 0;

  // Run the QLT algorithm with the values set by set_{rho,Q}. It is an error to
  // call this function from a parallel region.
  virtual void run() = 0;

  // Get a cell's tracer mass Qm after the QLT algorithm has run.
  virtual Real get_Qm(const Int& lclcellidx, const Int& tracer_idx) = 0;
};
} // namespace cedr

#endif
