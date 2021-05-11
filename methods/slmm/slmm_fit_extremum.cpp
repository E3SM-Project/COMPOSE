#include "slmm_fit_extremum.hpp"
#include "slmm_nla.hpp"
#include "slmm_gll.hpp"
#include "siqk_geometry.hpp"

namespace slmm {

class Projection {
  const Int m_, n_;
  std::vector<Real> Mtgt_, Mmix_;

public:
  // Mtgt is mxm and Mmix is mxn.
  Projection (const Int& m, const Int& n) : m_(m), n_(n) {
    Mtgt_.resize(m*m);
    Mmix_.resize(m*n);
  }

  Real* get_Mtgt () { return Mtgt_.data(); }
  Real* get_Mmix () { return Mmix_.data(); }

  void factor () {
#ifndef NDEBUG
    const int info =
#endif
      dpotrf('L', m_, Mtgt_.data(), m_);
    assert(info == 0);
  }

  void solve (const Real* const src, Real* const tgt) const {
    // tgt = Mmix src
    dgemm('n', 'n', m_, 1, n_, 1, Mmix_.data(), m_, src, n_, 0, tgt, m_);
    // tgt = Mtgt \ tgt
#ifndef NDEBUG
    const int info =
      dpotrs('L', m_, 1, Mtgt_.data(), m_, tgt, m_);
#endif
    assert(info == 0);
  }
};

/*
  Fit np=3 to the user's np function. This is equivalent to fitting
      m(x,y) := a22 x^2 y^2
              + b21 x^2 y + b12 x y^2
              + c20 x^2 + c11 x y + c02 y^2
              + d10 x + d01 y
              + e,
  for which
      m_x = 2 (a22 y^2 + b21 y + c20) x + (b12 y^2 + c11 y + d10)
      m_y = 2 (a22 x^2 + b12 x + c02) y + (b21 x^2 + c11 x + d01)
      m_xx = 2 (a22 y^2 + b21 y + c20)
      m_xy = 2 (2 a22 x + b12) y + (2 b21 x + c11)
           = 4 a22 x y + 2 (b21 x + b12 y) + c11
      m_yy = 2 (a22 x^2 + b12 x + c02).

  coef = [a22 b21 b12 c20 c11 c02 d10 d01 e]
           0   1   2   3   4   5   6   7  8
 */
class FitExtremum::Impl {
  static const constexpr Real
    bound_1d = 1,
    bound_2d = 1,
    max_relerr_1d_np6 = 0.025,
    max_relerr_2d_np6 = max_relerr_1d_np6;

  Int np_, np2_;
  Real max_relerr_1d, max_relerr_2d;

  struct Diagnostics {
    Real relerrs[5];
  };

  // L2 projection.
  std::shared_ptr<Projection> proj_2d_, proj_1d_;

  static Real eval_2d_basis (const Int& k, const Real& x, const Real& y) {
    switch (k) {
    case 0: return x*x*y*y;
    case 1: return x*x*y;
    case 2: return x*y*y;
    case 3: return x*x;
    case 4: return x*y;
    case 5: return y*y;
    case 6: return x;
    case 7: return y;
    case 8: return 1;
    default: assert(0); return 0;
    }
  }

  static Real eval_1d_basis (const Int& k, const Real& x) {
    switch (k) {
    case 0: return x*x;
    case 1: return x;
    case 2: return 1;
    default: assert(0); return 0;
    }
  }

  void init_fit () {
    using PlaneGeo = siqk::PlaneGeometry;

    const Int np = np_;
    max_relerr_1d = (np/6.0)*max_relerr_1d_np6;
    max_relerr_2d = (np*np/36.0)*max_relerr_2d_np6;

    GLL gll;
    const Real* gx, * gw;
    gll.get_coef(np, gx, gw);

    proj_2d_ = std::make_shared<Projection>(9, np*np);
    Real* Mtgt = proj_2d_->get_Mtgt();
    Real* Mmix = proj_2d_->get_Mmix();
    {
      siqk::TriangleQuadrature tq;
      siqk::RawConstVec3s tq_bary;
      siqk::RawConstArray tq_w;
      tq.get_coef(20, tq_bary, tq_w);
      const Int nq = len(tq_w);
      const Real* const ps = siqk::sqr::get_ref_vertices();
      for (Int r = 0; r < 9; ++r) {
        for (Int k = 1; k <= 2; ++k) {
          for (Int q = 0; q < nq; ++q) {
            Real rx[2];
            PlaneGeo::bary2coord(ps, ps+2*k, ps+2*(k+1), slice(tq_bary, q), rx);
            Real gj[GLL::np_max], gi[GLL::np_max]; {
              gll.eval(np, rx[0], gj);
              gll.eval(np, rx[1], gi);
            }
            const Real f = eval_2d_basis(r, rx[0], rx[1]);
            for (Int i = 0, c = 0; i < np; ++i)
              for (Int j = 0; j < np; ++j, ++c)
                Mmix[9*c + r] += tq_w[q]*f*gj[j]*gi[i];
            for (Int c = 0; c <= r; ++c)
              Mtgt[9*c + r] += tq_w[q]*f*eval_2d_basis(c, rx[0], rx[1]);
          }
        }
      }
    }
    proj_2d_->factor();

    proj_1d_ = std::make_shared<Projection>(3, np);
    Mtgt = proj_1d_->get_Mtgt();
    Mmix = proj_1d_->get_Mmix();
    {
      const Int nq = 12;
      const Real* qx, * qwt;
      gll.get_coef(nq, qx, qwt);
      for (Int r = 0; r < 3; ++r) {
        for (Int q = 0; q < nq; ++q) {
          Real g[GLL::np_max];
          gll.eval(nq, qx[q], g);
          const Real f = eval_1d_basis(r, qx[q]);
          for (Int c = 0; c < np; ++c)
            Mmix[3*c + r] += qwt[q]*f*g[c];
          for (Int c = 0; c <= r; ++c)
            Mtgt[3*c + r] += qwt[q]*f*eval_1d_basis(c, qx[q]);
        }
      }
    }
    proj_1d_->factor();
  }

  Real fit_2d (const Real* const y_gll, Real coef[9],
               Diagnostics& d) const {
    proj_2d_->solve(y_gll, coef);
    Real relerr; {
      GLL gll;
      const Real* gx, * gw;
      gll.get_coef(np_, gx, gw);
      Real num = 0, den = 0;
      for (Int i = 0; i < np_; ++i)
        for (Int j = 0; j < np_; ++j) {
          const Real x = gx[j], y = gx[i];
          const Real f = eval(coef, x, y) - coef[2];
          Real g = y_gll[np_*i+j] - coef[2];
          num += gw[i]*(f - g)*(f - g);
          den += gw[i]*g*g;
        }
      relerr = std::sqrt(num/den);
      d.relerrs[4] = relerr;
    }
    return relerr;
  }

  Real fit_1d (const Real* const y_gll, const Int dir, Real coef[3],
               Diagnostics& diag) const {
    Real d[GLL::np_max];
    switch (dir) {
    case 0: for (Int i = 0; i < np_; ++i) d[i] = y_gll[np_-1 + np_*i]; break;
    case 1: for (Int i = 0; i < np_; ++i) d[i] = y_gll[np_*(np_-1) + i]; break;
    case 2: for (Int i = 0; i < np_; ++i) d[i] = y_gll[np_*i]; break;
    case 3: for (Int i = 0; i < np_; ++i) d[i] = y_gll[i]; break;
    default: assert(0);
    }
    proj_1d_->solve(d, coef);
    Real relerr; {
      GLL gll;
      const Real* gx, * gw;
      gll.get_coef(np_, gx, gw);
      Real num = 0, den = 0;
      for (Int i = 0; i < np_; ++i) {
        const Real x = gx[i];
        const Real f = (coef[0]*x + coef[1])*x;
        Real g = d[i] - coef[2];
        num += gw[i]*(f - g)*(f - g);
        den += gw[i]*g*g;
      }
      relerr = std::sqrt(num/den);
      diag.relerrs[dir] = relerr;
    }
    return relerr;
  }

  void init (const Int np) {
    np_ = np;
    np2_ = np*np;
    init_fit();
  }

  void calc_gradient (const Real c[9], const Real& x, const Real& y,
                      Real grad[2]) const {
    grad[0] = 2*(c[0]*y*y + c[1]*y + c[3])*x + (c[2]*y*y + c[4]*y + c[6]);
    grad[1] = 2*(c[0]*x*x + c[2]*x + c[5])*y + (c[1]*x*x + c[4]*x + c[7]);
  }

  // H = {a b c} = [a b; b c]
  void calc_hessian (const Real c[9], const Real& x, const Real& y,
                     Real H[3]) const {
    H[0] = 2*(c[0]*y*y + c[1]*y + c[3]);
    H[1] = 4*c[0]*x*y + 2*(c[1]*x + c[2]*y) + c[4];
    H[2] = 2*(c[0]*x*x + c[2]*x + c[5]);
  }

  void calc_eigs (const Real H[3], Real eig[2]) const {
    const Real& a = H[0];
    const Real& b = H[1];
    const Real& c = H[2];
    const Real disc_sqrt = std::sqrt((a - c)*(a - c) + 4*b*b);
    eig[0] = 0.5*((a + c) + disc_sqrt);
    eig[1] = 0.5*((a + c) - disc_sqrt);
  }

  Real eval (const Real c[9], const Real& x, const Real& y) const {
    const Real x2 = x*x, y2 = y*y;
    return (c[0]*x2*y2 +
            c[1]*x2*y + c[2]*x*y2 +
            c[3]*x2 + c[4]*x*y + c[5]*y2 +
            c[6]*x + c[7]*y + c[8]);
  }

  bool calc_quadratic_extremum (const Real* y_gll, const Int dir,
                                Real& value, Diagnostics& d) const {
    Real coef[3];
    const Real re = fit_1d(y_gll, dir, coef, d);
    if (re > max_relerr_1d) return false;
    const Real a = coef[0], b = coef[1], c = coef[2];
    if (a == 0) return false;
    const Real x = -b/(2*a);
    if (x < -bound_1d || x > bound_1d) return false;
    value = (a*x + b)*x + c;
    return true;
  }

  bool solve_on_edges (const Real* y_gll, Real& min, Real& max,
                       Diagnostics& d) const {
    bool use = false;
    Real value;
    const auto update_extrema = [&] (const Real& value) {
      if (use) {
        min = std::min(min, value);
        max = std::max(max, value);
      } else {
        min = max = value;
        use = true;
      }
    };
    for (Int dir = 0; dir < 4; ++dir)
      if (calc_quadratic_extremum(y_gll, dir, value, d))
        update_extrema(value);
    return use;
  }

  bool is_in_bounds (const Real& x, const Real& y) const {
    return x >= -bound_2d && x <= bound_2d && y >= -bound_2d && y <= bound_2d;
  }

  bool solve (const Real c[9], Real& x, Real& y) const {
    const Int max_it = 5;
    x = y = 0;
    for (Int it = 0; it < max_it; ++it) {
      Real g[2], H[3];
      calc_gradient(c, x, y, g);
      calc_hessian(c, x, y, H);
      const Real det = H[0]*H[2] - H[1]*H[1];
      if (det == 0) { x = y = 2; return false; }
      assert(det != 0);
      const Real dx = (-H[2]*g[0] + H[1]*g[1])/det;
      const Real dy = ( H[1]*g[0] - H[0]*g[1])/det;
      x += dx;
      y += dy;
    }
    return is_in_bounds(x, y);
  }

  bool is_definite (const Real coef[9], Real& sign) const {
    sign = 0;
    bool definite = true;
    for (const Real x : {-1, 0, 1}) {
      for (const Real y : {-1, 0, 1}) {
        Real H[3], eig[2];
        calc_hessian(coef, x, y, H);
        calc_eigs(H, eig);
        for (int i = 0; i < 2; ++i)
          if (std::abs(eig[i]) > 1e2*std::numeric_limits<Real>::epsilon()) {
            const Real s = eig[i] > 0 ? 1 : -1;
            if (sign == 0) sign = s;
            if (sign != 0 && s != sign) {
              definite = false;
              break;
            }
          }
        if ( ! definite) break;
      }
      if ( ! definite) break;
    }
    return sign != 0 && definite;
  }

public:
  Impl (const Int np) {
    init(np);
  }

  bool calc_extrema (const Real* y_gll, Real coef[9], Real& min, Real& max) const {
    Diagnostics d;
    const bool use_1d = solve_on_edges(y_gll, min, max, d);
    Real sign;
    bool use_2d = (fit_2d(y_gll, coef, d) <= max_relerr_2d &&
                   (true || is_definite(coef, sign)));
    {
      bool all_relerr_bounded = true;
      for (Int i = 0; i < 4; ++i)
        if (d.relerrs[i] > max_relerr_1d)
          all_relerr_bounded = false;
      if (d.relerrs[4] > max_relerr_2d) all_relerr_bounded = false;
      if ( ! all_relerr_bounded) return false;
    }
    if ( ! use_2d) return use_1d;
    Real xi, yi;
    use_2d = solve(coef, xi, yi);
    if (use_2d) {
      const Real iv = eval(coef, xi, yi);
      if (use_1d) {
        min = std::min(min, iv);
        max = std::max(max, iv);
      } else {
        min = max = iv;
      }
    }
    return use_1d || use_2d;
  }
};

FitExtremum::FitExtremum (const Int np) {
  impl_ = std::make_shared<Impl>(np);
}

void FitExtremum::calc (const Real* const y_gll, Real& min, Real& max, bool& use,
                        Real* user_coef) {
  Real coef[9];
  use = impl_->calc_extrema(y_gll, coef, min, max);
  if (user_coef)
    for (Int i = 0; i < 9; ++i)
      user_coef[i] = coef[i];
}

} // namespace slmm
