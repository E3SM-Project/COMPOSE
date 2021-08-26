#ifndef INCLUDE_ISLET_MAXEIGCOMP_HPP
#define INCLUDE_ISLET_MAXEIGCOMP_HPP

#include "islet_types.hpp"
#include "islet_util.hpp"
#include "islet_interpmethod.hpp"
#include "islet_npx.hpp"

using islet::Array;

#include <memory>

class MaxEigComputer {
  static constexpr int max_nthread = 136;
  
  struct Workspace {
    Array<Real> A, wr, wi, work;
  };
  std::vector<Workspace> wss;
  Int max_ne;
  bool threaded, bloch;

public:
  MaxEigComputer (const bool ithreaded = true, const bool ibloch = true)
    : max_ne(0), threaded(ithreaded), bloch(ibloch)
  { setup_workspace(); }

  bool is_threaded () const { return threaded; }

  void setup_workspace(const Int max_ne_ = 64);

  Real run(const Int& ne_max, const Int& ndx_max, const Real& maxeigampm1,
           const bool quiet, const InterpMethod& im);

  Real run (const Int& np, const Int& ne_max, const Int& ndx_max,
            const Real& maxeigampm1, const bool quiet,
            UserInterpMethod* uim) {
    std::shared_ptr<UserInterpMethod> suim(uim, [] (UserInterpMethod*) {});
    return run(ne_max, ndx_max, maxeigampm1, quiet, InterpMethod(suim));
  }
  Real run (const Int& np, const Int& ne_max, const Int& ndx_max,
            const Real& maxeigampm1, const bool quiet,
            const UserInterpMethod::Ptr& uim) {
    return run(ne_max, ndx_max, maxeigampm1, quiet, InterpMethod(uim));
  }

  struct Analysis {
    static constexpr Real condv_switch = 1e2;
    Real max_eig_amp, max_condv, max_defect_ub;
  };

  Analysis calc_max_vals(const Int& nmu, const Int& ndx,
                         const InterpMethod& im);

  Analysis calc_max_vals (const Int& nmu, const Int& ndx, const Int& np,
                          UserInterpMethod* uim) {
    std::shared_ptr<UserInterpMethod> suim(uim, [] (UserInterpMethod*) {});
    return calc_max_vals(nmu, ndx, InterpMethod(suim));
  }

  // dx is fraction of an element, so [0,1], *not* reference-element cordinate.
  void compute(const InterpMethod& im, const Real& dx, const Int& ne,
               Real* max_amp_out, Real* max_condv = nullptr,
               Real* max_defect_ub = nullptr,
               Complex* lam = nullptr, Complex* V = nullptr);

  void compute(UserInterpMethod* uim, const Real& dx, const Int& ne,
               Real* max_amp_out, Real* max_condv = nullptr,
               Real* max_defect_ub = nullptr,
               Complex* lam = nullptr, Complex* V = nullptr) {
    std::shared_ptr<UserInterpMethod> suim(uim, [] (UserInterpMethod*) {});
    compute(InterpMethod(suim), dx, ne, max_amp_out, max_condv,
            max_defect_ub, lam, V);
  }

  static int unittest();
};

#endif
