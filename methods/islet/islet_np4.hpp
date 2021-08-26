#ifndef INCLUDE_ISLET_NP4_HPP
#define INCLUDE_ISLET_NP4_HPP

#include "islet_types.hpp"
#include "islet_isl.hpp"
#include "islet_interpmethod.hpp"

struct Np4InterpMethod : public UserInterpMethod {
  typedef std::shared_ptr<Np4InterpMethod> Ptr;

  Np4InterpMethod(Real c0, Real c1, Real c2);
  void reset_c(Real c0, Real c1, Real c2);
  Real eval_a(const Real& x) const;
  void eval(const Real& x, Real* const y);
  const Real* get_xnodes () const override { return islet::get_x_gll(4); }
  Int get_np () const override { return 4; }

private:
  Real c[3];
};

#endif
