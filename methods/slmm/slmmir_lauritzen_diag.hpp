#ifndef INCLUDE_SLMMIR_LAURITZEN_DIAG_HPP
#define INCLUDE_SLMMIR_LAURITZEN_DIAG_HPP

#include <memory>

#include "slmm_gallery.hpp"
#include "slmmir_remap_data.hpp"
#include "slmmir_d2c.hpp"

struct LauritzenDiag {
  typedef std::shared_ptr<LauritzenDiag> Ptr;

  LauritzenDiag(const Int nsteps_per_12days, const Int len,
                const std::vector<gallery::InitialCondition::Shape>& ics,
                const Real* tracer_data, const Real* dA,
                const bool expensive_io);

  bool run(const Int step, const Real* tracer_data, const Real* dA,
           const D2Cer& d2cer, const std::string& out_prefix);
  void run_me(const Int len, const Real* tracer_data, const Real* dA,
              const bool first);
  void print();

private:
  class Impl;
  std::shared_ptr<Impl> p;
};

#endif
