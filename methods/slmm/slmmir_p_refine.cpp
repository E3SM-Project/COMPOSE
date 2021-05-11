#include "slmmir_p_refine.hpp"
#include "slmmir_util.hpp"

Transferer2D
::Transferer2D (Int np_from_, const Basis::Ptr& b_from_, Int np_to_, const Basis::Ptr& b_to_) {
  np_from = np_from_; np_to = np_to_;
  np_from2 = square(np_from); np_to2 = square(np_to);
  b_from = b_from_; b_to = b_to_;
}

MatrixBasedTransferer2D
::MatrixBasedTransferer2D (Int np_from, const Basis::Ptr& b_from, Int np_to, const Basis::Ptr& b_to)
  : Transferer2D(np_from, b_from, np_to, b_to)
{}

void MatrixBasedTransferer2D::apply (const Real* src, Real* dst) const {
  for (Int i = 0; i < np_to2; ++i) {
    const Real* opi = op.data() + np_from2*i;
    Real accum = 0;
    for (Int j = 0; j < np_from2; ++j)
      accum += opi[j]*src[j];
    dst[i] = accum;
  }
}

Interpolator2D
::Interpolator2D (Int np_from, const Basis::Ptr& b_from, Int np_to, const Basis::Ptr& b_to)
  : MatrixBasedTransferer2D(np_from, b_from, np_to, b_to)
{
  init();
}

void Interpolator2D::init () {
  op.optclear_and_resize(np_from2*np_to2);
  const Real* x_gll_to;
  b_to->get_x(np_to, x_gll_to);
  Real fex[GLL::np_max], fey[GLL::np_max];
  for (Int tj = 0; tj < np_to; ++tj) {
    if (tj == 0 || tj == np_to-1) {
      for (Int i = 0; i < np_from; ++i) fey[i] = 0;
      fey[tj == 0 ? 0 : np_from-1] = 1;
    } else {
      b_from->eval(np_from, x_gll_to[tj], fey);
    }
    for (Int ti = 0; ti < np_to; ++ti) {
      if (ti == 0 || ti == np_to-1) {
        for (Int i = 0; i < np_from; ++i) fex[i] = 0;
        fex[ti == 0 ? 0 : np_from-1] = 1;
      } else {
        b_from->eval(np_from, x_gll_to[ti], fex);
      }
      for (Int fj = 0; fj < np_from; ++fj)
        for (Int fi = 0; fi < np_from; ++fi)
          op[(tj*np_to + ti)*np_from2 + fj*np_from + fi] = fex[fi]*fey[fj];
    }
  }
}

void MeshInterpolator
::apply (const Mesh& m_coarse, const Mesh& m_fine,
         const AVec3s& p_coarse,
         AVec3s& p_fine) const {
  const Int ne = nslices(m_coarse.dglln2cglln) / square(m_coarse.np);
  assert(ne == nslices(m_fine.dglln2cglln) / square(m_fine.np));
  const Int np_from = i2d.get_np_from(), np_to = i2d.get_np_to();
# pragma omp parallel for
  for (Int ie = 0; ie < ne; ++ie) {
    // Interpolate coarse points to fine within the element using the standard
    // GLL interpolant.
    const Int* const cc = &m_coarse.dglln2cglln(ie*square(np_from));
    const Int* const cf = &m_fine.dglln2cglln(ie*square(np_to));
    Real from[GLL::np_max*GLL::np_max], to[GLL::np_max*GLL::np_max];
    for (Int d = 0; d < 3; ++d) {
      for (Int cj = 0; cj < np_from; ++cj)
        for (Int ci = 0; ci < np_from; ++ci) {
          const auto v = slice(p_coarse, cc[cj*np_from + ci]);
          from[cj*np_from + ci] = v[d];
        }
      i2d.apply(from, to);
      for (Int fj = 0; fj < np_to; ++fj)
        for (Int fi = 0; fi < np_to; ++fi) {
          auto wv = slice(p_fine, cf[fj*np_to + fi]);
          wv[d] = to[fj*np_to + fi];
        }
    }
  }
}

/*
  B x + B B'lam = B y
  B B'lam = B y - u
  lam = (B B')"(B y - u)
  x = y - B'lam
    = y - B'(B B')"(B y - u)
  B' = op = Q R
  B B' = R'Q'Q R = R'R
  B'(B B')"B = Q R (R'R)"R'Q' = Q Q'
  B'(B B')" = Q R (R'R)" = Q R'"
  x = y - Q Q'y + Q R'\u
    = y + Q (R'\u - Q'y)
*/
void Interpolator2D
::find_x_nearest_y_interpolating_u (const Real* y, const Real* u, Real* x) {
  assert(0); // We're not using this function right now.
  if ( ! qr_fac) {
#   pragma omp barrier
#   pragma omp master
    {
      qr_fac = std::make_shared<QrFac>(np_from2, np_to2);
      std::copy(op.data(), op.data() + np_from2*np_to2, qr_fac->A());
      qr_fac->factor();
    }
#   pragma omp barrier
  }
  // wrk = R'\u
  Real wrk[GLL::np_max*GLL::np_max];
  std::copy(u, u + np_to2, wrk);
  qr_fac->solve_Rt(wrk);
  // x = Q'y
  qr_fac->apply_Qt(y, x);
  // wrk = R'\u - Q'y
  for (Int i = 0; i < np_to2; ++i) wrk[i] -= x[i];
  // x = Q (R'\u - Q'y)
  qr_fac->apply_Q(wrk, x);
  // and done
  for (Int i = 0; i < np_from2; ++i) x[i] += y[i];
}

void calc_pref_gll_quantities (const Mesh& mc, const RemapData& rdc, const Mesh& mf,
                               RemapData& rdf) {
  Interpolator2D c2f(mc.np, mc.basis, mf.np, mf.basis);
  const Int npc2 = square(mc.np), npf2 = square(mf.np), nelem = nslices(mc.geo_c2n);
# pragma omp parallel for
  for (Int ic = 0; ic < nelem; ++ic)
    c2f.apply(&rdc.Jt_[ic*npc2], &rdf.Jt_[ic*npf2]);
  calc_gll_basis_function_integrals(mf.np, *mf.basis, mf.geo_c2n, mf.geo_p,
                                    rdf.Jt_.data(), rdf.dgbfi_gll_.data());
}
