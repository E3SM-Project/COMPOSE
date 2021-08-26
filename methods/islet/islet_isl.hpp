#ifndef INCLUDE_ISLET_ISL_HPP
#define INCLUDE_ISLET_ISL_HPP

#include "islet_types.hpp"
#include "islet_interpmethod.hpp"

#include <memory>

namespace islet {
const Real* get_x_gll(const Int np);
const Real* get_w_gll(const Int np);

template <typename Scalar>
void eval_lagrange_poly (const Scalar* x_gll, const Int& np, const Scalar& x,
                         Scalar* const y) {
  for (int i = 0; i < np; ++i) {
    Scalar f = 1;
    for (int j = 0; j < np; ++j)
      f *= (i == j) ?
        1 :
        (x - x_gll[j]) / (x_gll[i] - x_gll[j]);
    y[i] = f;
  }
}

struct Operator {
  typedef std::shared_ptr<Operator> Ptr;
  typedef std::shared_ptr<const Operator> ConstPtr;

  virtual void eval(const Int& np, const Real& x, Real* const v) const = 0;
  virtual const Real* get_xnodes (const Int& np) const { return get_x_gll(np); }
  virtual std::string get_basis_string (const Int& np) const { return ""; }

  enum Method { gll_natural = 0, gll_offset_nodal_subset, xnodal, gll_best,
                uniform_offset_nodal_subset };
  static ConstPtr create(Method m);
};

struct OperatorInterpMethod : public UserInterpMethod {
  typedef std::shared_ptr<OperatorInterpMethod> Ptr;
  OperatorInterpMethod (const Int np_, const Operator::ConstPtr& op_) : np(np_), op(op_) {}
  void eval (const Real& x, Real* const v) override { op->eval(np, x, v); }
  const Real* get_xnodes () const override { return op->get_xnodes(np); }
  Int get_np () const override { return np; }
private:
  Int np;
  Operator::ConstPtr op;
};

Int unittest_eval();

} // namespace islet

#endif
