#ifndef INCLUDE_ISLET_XNODES_METRICS_HPP
#define INCLUDE_ISLET_XNODES_METRICS_HPP

#include "islet_types.hpp"
#include "islet_nodalbasis.hpp"

#include <vector>
#include <sstream>

// l1 only
Real calc_xnodes_metric(const Nodes& nodes, const Real* const xnodes);
// l1, l2, linf
void calc_xnodes_metrics(const Nodes& nodes, const Real* const xnodes, Real* metrics);

void calc_lebesgue_consts(const Nodes& nodes, const Real* const xnodes, Real* metrics);

void calc_weights(const Nodes& nodes, const Real* const xnode, Real* const wt);

struct MetricsTracker {
  typedef std::shared_ptr<MetricsTracker> Ptr;

  MetricsTracker(Int np, bool very_strict = false);

  void set_pum_max(Real pum); // optional; default is 1
  Real get_pum_max () const { return pum_max; }

  // Min pum seen so far. If none, return 1.
  Real get_pum_min () const { return pum_min; }

  // Compute metrics. Return whether these are provisionally acceptable.
  bool acceptable_metrics(const Nodes& nodes, const Real* xnodes,
                          const Real* metrics) const;
  // pum needs to be <= this value to update.
  Real pum_to_accept(const Nodes& nodes, const Real* xnodes,
                     const Real* metrics) const;
  // Would update based on metrics and pum?
  bool would_update(const Real* metrics, const Real& pum) const;
  // Do the update.
  void update(const Real* metrics, const Real& pum);

  void get_metrics(Real pum, Real* metrics) const;

  bool write(std::ofstream& os);
  bool read(std::ifstream& os);

private:
  static const Int nmet = 3, nbin = 30;
  Real best_metrics[nmet*nbin]; // l1, l2, linf
  Real pum_bins[nbin+1];
  Real pum_max, pum_min;
};

#endif
