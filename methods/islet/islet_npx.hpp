#ifndef INCLUDE_ISLET_NPX_HPP
#define INCLUDE_ISLET_NPX_HPP

#include "islet_util.hpp"
#include "islet_types.hpp"
#include "islet_tables.hpp"
#include "islet_interpmethod.hpp"

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

template <typename Scalar=Real>
struct npx {
  static void eval (const Int np, const Scalar& x, Scalar* const v) {
    eval_lagrange_poly(islet::get_x_gll(np), np, x, v);
  }
};

template <typename Scalar=Real>
struct npxstab : public npx<Scalar> {
  template <int np, int nreg>
  static void eval (const std::array<int, nreg>& order,
                    const std::array<int, nreg>& os,
                    const Scalar& x, Scalar* const v) {
    if (x > 0) {
      eval<np, nreg>(order, os, -x, v);
      for (int i = 0; i < np/2; ++i)
        std::swap(v[i], v[np-i-1]);
      return;
    }
    const auto x_gll = islet::get_x_gll(np);
    bool done = false;
    for (Int i = 0; i < nreg; ++i)
      if (x < x_gll[i+1]) {
        std::fill(v, v + np, 0);
        eval_lagrange_poly(x_gll + os[i], order[i], x, v + os[i]);
        done = true;
        break;
      }
    if ( ! done)
      eval_lagrange_poly(x_gll, np, x, v);
  }

  static void eval(const Int& np, const Scalar& x, Scalar* const v);

  static int ooa_vs_np (const int np) {
    if (np == 5) return 2;
    return np - 1 - ((np-1)/3);
  }

  template <int np, int nreg, int alphanp>
  static void eval (const std::array<int, nreg>& subnp,
                    const std::array<int, nreg>& os,
                    const std::array<Scalar, nreg*alphanp>& alphac,
                    const Scalar& x, Scalar* const v) {
    if (x > 0) {
      eval<np, nreg, alphanp>(subnp, os, alphac, -x, v);
      for (int i = 0; i < np/2; ++i)
        std::swap(v[i], v[np-i-1]);
      return;
    }
    const auto x_gll = islet::get_x_gll(np);
    bool done = false;
    for (int i = 0; i < nreg; ++i)
      if (x < x_gll[i+1]) {
        eval_lagrange_poly(x_gll, np, x, v);
        if (subnp[i] < np) {
          Real w[12] = {0};
          eval_lagrange_poly(x_gll + os[i], subnp[i], x, w + os[i]);
          Real alpha = 0;
          if (alphanp == 1)
            alpha = alphac[i];
          else {
            assert(alphanp <= 3);
            const auto alpha_r_gll = islet::get_x_gll(alphanp);
            const auto r = (x - x_gll[i]) / (x_gll[i+1] - x_gll[i]);
            Real a[3];
            eval_lagrange_poly(alpha_r_gll, alphanp, r, a);
            for (int j = 0; j < alphanp; ++j)
              alpha += alphac[alphanp*i + j] * a[j];
          }
          for (int j = 0; j < np; ++j)
            v[j] = alpha*v[j] + (1 - alpha)*w[j];
        }
        done = true;
        break;
      }
    if ( ! done)
      eval_lagrange_poly(x_gll, np, x, v);
  }
};

struct InterpMethod {
  typedef std::shared_ptr<InterpMethod> Ptr;

  enum Type { notype, npx, npxstab, user };
  Int np;
  Type type;
  std::shared_ptr<UserInterpMethod> uim;

  static Type convert (const std::string& s) {
    if (s == "npx") return npx;
    if (s == "npxstab") return npxstab;
    if (s == "user") return user;
    throw std::runtime_error(std::string("Not an InterpMethod::Type: ") + s);
  }

  static std::string convert (const Type& t) {
    if (t == npx) return "npx";
    if (t == npxstab) return "npxstab";
    if (t == user) return "user";
    throw std::runtime_error("Not an InterpMethod::Type.");
  }

  InterpMethod () : np(-1), type(notype) {}
  InterpMethod (Int inp, Type itype) : np(inp), type(itype) {}
  InterpMethod (const std::shared_ptr<UserInterpMethod>& iuim)
    : np(iuim->get_np()), type(user), uim(iuim) {}

  const Real* get_xnodes() const {
    return uim ? uim->get_xnodes() : islet::get_x_gll(np);
  };
};

void op_eval(const InterpMethod& im, const Real a_src, Real* v);

#endif
