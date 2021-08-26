#ifndef INCLUDE_ISLET_INTERPMETHOD_HPP
#define INCLUDE_ISLET_INTERPMETHOD_HPP

#include <memory>

struct UserInterpMethod {
  typedef std::shared_ptr<UserInterpMethod> Ptr;
  virtual ~UserInterpMethod () {}
  virtual void eval(const Real& x, Real* const v) = 0;
  virtual const Real* get_xnodes() const = 0;
  virtual Int get_np() const = 0;
};

#endif
