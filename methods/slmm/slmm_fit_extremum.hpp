#ifndef INCLUDE_SLMM_FIT_EXTREMUM_HPP
#define INCLUDE_SLMM_FIT_EXTREMUM_HPP

#include <memory>

#include "slmm_defs.hpp"

namespace slmm {

struct FitExtremum {
  FitExtremum(const Int np);

  void calc(const Real* const y_gll, Real& min, Real& max, bool& use,
            Real* coef = nullptr);

private:
  class Impl;
  std::shared_ptr<Impl> impl_;
};

} // namespace slmm

#endif
