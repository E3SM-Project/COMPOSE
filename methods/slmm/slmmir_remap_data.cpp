#include "slmmir_remap_data.hpp"

#include "slmmir_util.hpp"
#include "slmm_gll.hpp"

void FullMassMatrix::init (const int nelem, const int np) {
  np_ = np;
  Array<int> rowptr(nelem + 1), colidx(nelem);
  for (int i = 0; i < nelem; ++i) {
    rowptr[i] = i;
    colidx[i] = i;
  }
  rowptr[nelem] = nelem;
  m_.init(nelem, nelem, np2(), np2(), rowptr.data(), colidx.data());
  m_.zero();
  assert(m_.m() == np2() && m_.n() == np2());
  assert(m_.M() == m_.N() && m_.M() == nelem);
  assert(m_.blockrow(0) + np4() == m_.blockrow(1));
}

const double* FullMassMatrix::block (const int i) const {
  assert(m_.blockrow(i) - m_.blockrow(0) == i*np4());
  return m_.blockrow(i);
}
double* FullMassMatrix::block (const int i) {
  return const_cast<double*>(const_cast<const MT&>(m_).blockrow(i));
}

void FullMassMatrix::factor () {
# pragma omp parallel for
  for (int i = 0; i < m_.M(); ++i)
    factor(i);
}

void FullMassMatrix::factor (const int i) {
  const int n = np2();
  double* const d = block(i);
  const char uplo = 'L';
  const int info = dpotrf(uplo, n, d, n);
  if (info != 0) {
    fprintf(stderr, "M() %d i %d info %d\n", m_.M(), i, info);
    fprintf(stderr, "a = [");
    for (int c = 0; c < n; ++c) {
      for (int r = 0; r < n; ++r)
        fprintf(stderr, " %1.15e", d[n*c + r]);
      fprintf(stderr, ";");
    }
    fprintf(stderr, "];\n");
  }
  assert(info == 0);
}

void FullMassMatrix::solve (const int elem, double* const bx, const int nrhs,
                            const int ldbx) const {
  const int n = np2();
  const double* const d = block(elem);
#ifndef NDEBUG
  const int info =
#endif
    dpotrs('L', n, nrhs, d, n, bx, ldbx);
  assert(info == 0);
}

void FullMassMatrix::fwdsub (const int elem, double* const bx, const int nrhs,
                             const int ldbx) const {
  const int n = np2();
  const double* const d = block(elem);
  dtrsm('L', 'L', 'N', 'N', n, nrhs, 1, d, n, bx, ldbx);
}

void FullMassMatrix::bwdsub (const int elem, double* const bx, const int nrhs,
                             const int ldbx) const {
  const int n = np2();
  const double* const d = block(elem);
  dtrsm('L', 'L', 'T', 'N', n, nrhs, 1, d, n, bx, ldbx);
}

void FullMassMatrix
::solve_1eq_ls (const int elem, Real* x, const int nrhs, const int ldx,
                Real* c, const bool c_is_Linvc, const Real d) const {
  // The Lagrangian is
  //     Lag = 1/2 x' A'A x - b' A x - mu (c' x - d),
  // giving the system
  //     Lag_x  = A'A x - A' b - c mu = 0
  //     Lag_mu = c'x - d = 0,
  // or
  //     [A'A -c] [x ] = [A'b]
  //     [c'    ] [mu]   [d  ].
  // We have L = chol(A'A). Then
  //     [L          0] [L'  -inv(L) c     ] [x ] = [A'b]
  //     [c'inv(L')  I] [0    c'inv(L L') c] [mu]   [d  ].
  // Hence
  //     1. [L          0] [a1] = [A'b]
  //        [c'inv(L')  I] [a2]   [d  ]
  //     => 1a. Solve L s = c
  //        1b. Solve L a1 = A'b
  //        1c. Set a2 = d - s'a1.
  //     2. [L'  -inv(L) c     ] [x ] = [a1]
  //        [0    c'inv(L L') c] [mu]   [a2]
  //     => 2a. Set mu = a2 / (s's)
  //        2b. Solve L' x = a1 + mu s.
  const int n = np2();
  /* 1a */ if ( ! c_is_Linvc) fwdsub(elem, c, 1, n);
  /* 1b */ fwdsub(elem, x, nrhs, ldx);
  /* 1c */
  Real a2 = 0;
  for (int i = 0; i < n; ++i) a2 += c[i]*x[i];
  a2 = d - a2;
  /* 2a */
  Real s2 = 0;
  for (int i = 0; i < n; ++i) s2 += c[i]*c[i];
  const Real mu = a2/s2;
  /* 2b */
  for (int i = 0; i < n; ++i) x[i] += mu*c[i];
  bwdsub(elem, x, nrhs, ldx);
}

class CalcM {
  const AIdxs& c2n_;
  FullMassMatrix& fmm_;
  const Int np_, np2_;
  const bool facet_transport_;
  siqk::TriangleQuadrature tq_;
  siqk::RawConstVec3s tq_bary_;
  siqk::RawConstArray tq_w_;
  Int nq_;
  GLL gll_;

  // Set by init().
  const AVec3s p_;

public:
  CalcM (const AIdxs& c2n, const Int& np,
         const Int tq_order, FullMassMatrix& fmm, const Dmc::Enum dmc)
    : c2n_(c2n), fmm_(fmm), np_(np), np2_(square(np_)),
      facet_transport_(Dmc::is_facet(dmc))
  {
    tq_.get_coef(tq_order, tq_bary_, tq_w_);
    nq_ = len(tq_w_);
  }

  // Must call init before op().
  void init (const AVec3s& p) {
    *const_cast<AVec3s*>(&p_) = p;
  }

  // Calculate for cell ci.
  void operator() (const Int ci) {
    if (facet_transport_) {
      // This is a lot of wasted work. But keep this to match QOS case.
      const Real* const ps = siqk::sqr::get_ref_vertices();
      Real* const block = fmm_.block(ci);
      for (Int k = 1; k <= 2; ++k) {
        const Real jac_reftri2refsquare = PlaneGeo::calc_tri_jacobian(
          ps, ps+2*k, ps+2*(k+1));
        for (Int q = 0; q < nq_; ++q) {
          Real square_coord[2];
          PlaneGeo::bary2coord(ps, ps+2*k, ps+2*(k+1), slice(tq_bary_, q),
                               square_coord);
          Real gj[GLL::np_max], gi[GLL::np_max]; {
            gll_.eval(np_, square_coord[0], gi);
            gll_.eval(np_, square_coord[1], gj);
          }
          const Real d0 = 0.5 * tq_w_[q] * jac_reftri2refsquare;
          for (Int aj = 0, a_basis_idx = 0; aj < np_; ++aj) {
            const Real d1 = d0 * gj[aj];
            for (Int ai = 0; ai < np_; ++ai, ++a_basis_idx) {
              Real d2 = d1 * gi[ai];
              for (Int bj = 0, b_basis_idx = 0; bj < np_; ++bj) {
                const Real d3 = d2 * gj[bj];
                for (Int bi = 0; bi < np_; ++bi, ++b_basis_idx) {
                  if (b_basis_idx < a_basis_idx) continue;
                  const Real d = d3 * gi[bi];
                  block[np2_*a_basis_idx + b_basis_idx] += d;
                  if (a_basis_idx != b_basis_idx)
                    block[np2_*b_basis_idx + a_basis_idx] += d;
                }
              }
            }
          }
        }
      }
    } else {
      Real ps[12];
      copy_vertices(p_, c2n_, ci, ps);
      const auto cell = slice(c2n_, ci);
      Real* const block = fmm_.block(ci);
      for (Int k = 1; k <= 2; ++k) {
        for (Int q = 0; q < nq_; ++q) {
          Real sphere_coord[3];
          const Real jac_reftri2sphere = SphereGeo::calc_tri_jacobian(
            ps, ps+3*k, ps+3*(k+1), slice(tq_bary_, q), sphere_coord);
          Real gj[GLL::np_max], gi[GLL::np_max]; {
            Real a, b;
            siqk::sqr::calc_sphere_to_ref(p_, cell, sphere_coord, a, b);
            gll_.eval(np_, b, gj);
            gll_.eval(np_, a, gi);
          }
          Real d0 = 0.5 * tq_w_[q] * jac_reftri2sphere;
          for (Int aj = 0, a_basis_idx = 0; aj < np_; ++aj) {
            const Real d1 = d0 * gj[aj];
            for (Int ai = 0; ai < np_; ++ai, ++a_basis_idx) {
              Real d2 = d1 * gi[ai];
              for (Int bj = 0, b_basis_idx = 0; bj < np_; ++bj) {
                const Real d3 = d2 * gj[bj];
                for (Int bi = 0; bi < np_; ++bi, ++b_basis_idx) {
                  if (b_basis_idx < a_basis_idx) continue;
                  const Real d = d3 * gi[bi];
                  block[np2_*a_basis_idx + b_basis_idx] += d;
                  if (a_basis_idx != b_basis_idx)
                    block[np2_*b_basis_idx + a_basis_idx] += d;
                }
              }
            }
          }
        }
      }
    }
  }
};

static void calc_M_fwd (const Mesh& m, FullMassMatrix& fmm,
                        const Dmc::Enum dmc) {
  fmm.init(nslices(m.geo_c2n), m.np);
  CalcM calc_m(m.geo_c2n, m.np, m.tq_order, fmm, dmc);
  calc_m.init(m.geo_p);
# pragma omp parallel for
  for (Int ci = 0; ci < nslices(m.geo_c2n); ++ci)
    calc_m(ci);
  fmm.factor();
}

RemapData::RemapData (const Mesh& m, const Method::Enum method, const Dmc::Enum dmc) {
  ot_.init(a2ConstVec3s(m.geo_p), a2ConstIdxs(m.geo_c2n));
  calc_node_jacobians(m, m.geo_p, Jt_);
  calc_basis_function_integrals(m, m.geo_p, dgbfi_, cgbfi_);
  calc_gll_basis_function_integrals(m, GLL(), dgbfi_gll_);
  method_ = method;
  dmc_ = dmc;
  if ( ! Method::is_pisl(method)) calc_M_fwd(m, fmm_, dmc);
  init_dgbfi_mass(Dmc::use_homme_mass(dmc_) ? dgbfi_gll_ : dgbfi_);
}

void RemapData
::init_dgbfi_mass (const ARealArray& F) {
  *const_cast<ARealArray*>(&dgbfi_mass_) = F;

  if (Dmc::is_equality_constrained(dmc_)) {
    resize(Linv_dgbfi_mass_, F.size());
    if (Dmc::is_facet(dmc_)) {
      // In the facet case, J_ref^sphere needs to be removed because it is part
      // of the field ({rho,Q} J).
      const Int dnn = F.size();
#     pragma omp parallel for
      for (Int i = 0; i < dnn; ++i)
        Linv_dgbfi_mass_(i) = F(i) / Jt_(i);
    } else {
      copy(Linv_dgbfi_mass_, F);
    }
    const Int ncell = fmm_.get_M().M();
#   pragma omp parallel for
    for (Int ci = 0; ci < ncell; ++ci)
      fmm_.fwdsub(ci, &Linv_dgbfi_mass_(ci*fmm_.np2()), 1, F.size());
  }
}

void RemapData::solve_M_full (Real* bx, const int nrhs,
                              const int ldxb) const {
# pragma omp parallel for
  for (MT::Int br = 0; br < T_.M(); ++br)
    fmm_.solve(br, bx + br*fmm_.np2(), nrhs, ldxb);
}

void RemapData::apply_T_cell (const int ti, const Real* src, const int lds,
                              Real* tgt, const int ldt, const int nrhs,
                              const Real* FsmoFtm) const {
  const Int np2 = T_.n(), np4 = square(np2);
  const auto rowptr = T_.rowptr();
  const auto colidx = T_.colidx();
  const RemapData::MT::Scalar* const values = T_.blockrow(0);
  Real* const y = tgt + np2*ti;
  for (Int k = 0; k < np2; ++k) y[k] = 0;
  if (Dmc::is_facet(dmc_)) {
    for (Int j = rowptr[ti]; j < rowptr[ti+1]; ++j) {
      const Int si = colidx[j];
      const auto os = si*np2;
      // Source field from transport method's perspective must be (Q J), where J
      // is the Jacobian of ref cell to sphere. Bring J in here.
      const Real* const Jt = Jt_.data() + os;
      Real src_si[GLL::np_max*GLL::np_max];
      for (Int i = 0; i < np2; ++i)
        src_si[i] = src[os + i]*Jt[i];
      // If we're using IR, multiply the source by the density factory. To be
      // slightly more efficient, we could put J into FsmoFtm. For clarity, we
      // don't.
      if (Method::is_ir(method_)) {
        const Real* const FsmoFtm0 = FsmoFtm + os;
        for (Int i = 0; i < np2; ++i)
          src_si[i] *= FsmoFtm0[i];
      }
      // Apply the mixed mass matrix.
      const RemapData::MT::Scalar* const b = values + j*np4;
      if (np2 == 4)
        matvec4(b, src_si, y, nrhs);
      else
        dgemm('t', 'n', np2, nrhs, np2, 1, b, np2, src_si, lds, 1, y, ldt);
    }
  } else {
    for (Int j = rowptr[ti]; j < rowptr[ti+1]; ++j) {
      const Int si = colidx[j];
      const RemapData::MT::Scalar* const b = values + j*np4;
      Real srclcl[GLL::np_max*GLL::np_max];
      const Real* const src0 = src + si*np2;
      // Density factor.
      if (Method::is_ir(method_)) {
        const Real* const FsmoFtm0 = FsmoFtm + si*np2;
        for (Int i = 0; i < np2; ++i)
          srclcl[i] = src0[i] * FsmoFtm0[i];
      }
      // Apply the mixed mass matrix.
      if (np2 == 4)
        matvec4(b,
                Method::is_ir(method_) ? srclcl : src0,
                y, nrhs);
      else
        dgemm('t', 'n', np2, nrhs, np2, 1, b, np2,
              Method::is_ir(method_) ? srclcl : src0,
              lds, 1, y, ldt);
    }
  }
}

void RemapData
::remap_cell (const Int ti, const Real* src, const int lds, Real* tgt,
              const int ldt, const int nrhs, const Real* FsmoFtm) const {
  if ( ! Dmc::is_equality_constrained(dmc_)) {
    const Int np2 = T_.n();
    apply_T_cell(ti, src, lds, tgt, ldt, nrhs, FsmoFtm);
    fmm_.solve(ti, tgt + np2*ti, nrhs, ldt);
  } else {
    const MT::Int n = T_.n(), tos = n*ti;
    apply_T_cell(ti, src, lds, tgt, ldt, nrhs, FsmoFtm);
    const auto rowptr = T_.rowptr();
    const auto colidx = T_.colidx();
    Real mass_t = 0;
    for (Int k = rowptr[ti]; k < rowptr[ti+1]; ++k) {
      const Real* const p_s_ol_c = p_s_ol_.data() + k*n;
      const auto sos = colidx[k]*n;
      const Real* const src_c = src + sos;
      const Real* const dgbfi_c = dgbfi_mass_.data() + sos;
      for (Int i = 0; i < n; ++i)
        mass_t += p_s_ol_c[i] * dgbfi_c[i] * src_c[i];
    }
    fmm_.solve_1eq_ls(ti, tgt + tos, 1, ldt,
                      const_cast<Real*>(Linv_dgbfi_mass_.data()) + tos,
                      true, mass_t);
  }
  if (Dmc::is_facet(dmc_)) {
    const Int np2 = T_.n();
    const auto os = np2*ti;
    Real* const tgt_ti = tgt + os;
    const Real* const Jt = Jt_.data() + os;
    for (Int i = 0; i < np2; ++i)
      tgt_ti[i] /= Jt[i];
  }
}

void RemapData::remap (const Real* src, const int lds, Real* tgt,
                       const int ldt, const int nrhs, const Real* FsmoFtm) const {
  const auto ncell = T_.M();
  if (dmc_ == Dmc::glbl_eq_homme) {
    const Int np2 = T_.n();
    const Int N = ncell*np2;
    Real mass = 0, a2 = 0, s2 = 0;
#   pragma omp parallel
    {
#     pragma omp for reduction (+:mass)
      for (Int i = 0; i < N; ++i)
        mass += dgbfi_mass_[i] * src[i];
#     pragma omp for
      for (Int ti = 0; ti < ncell; ++ti) {
        apply_T_cell(ti, src, lds, tgt, ldt, nrhs, FsmoFtm);
        fmm_.fwdsub(ti, tgt + np2*ti, nrhs, ldt);
      }
#     pragma omp for reduction (+:a2)
      for (Int i = 0; i < N; ++i)
        a2 += Linv_dgbfi_mass_[i] * tgt[i];
#     pragma omp single
      a2 = mass - a2;
#     pragma omp for reduction (+:s2)
      for (Int i = 0; i < N; ++i)
        s2 += square(Linv_dgbfi_mass_[i]);  //todo cache
      const Real mu = a2/s2;
#     pragma omp for
      for (Int i = 0; i < N; ++i)
        tgt[i] += mu*Linv_dgbfi_mass_[i];
#     pragma omp for
      for (Int ti = 0; ti < ncell; ++ti)
        fmm_.bwdsub(ti, tgt + np2*ti, nrhs, ldt);
    }
  } else {
#   pragma omp parallel for
    for (Int ti = 0; ti < ncell; ++ti)
      remap_cell(ti, src, lds, tgt, ldt, nrhs, FsmoFtm);
  }
}

static void report (const std::string label, const Real* const x_t,
                    const Real* const x, const Int n) {
  Real me = 0, den = 0;
  for (Int i = 0; i < n; ++i) {
    me = std::max(me, std::abs(x[i] - x_t[i]));
    den = std::max(den, std::abs(x_t[i]));
  }
  printf("> RemapData %21s: %1.3e\n", label.c_str(), me/den);
}

void RemapData::compare_MT () const {
  Real diag_num = 0, diag_den = 0;
  const auto& M = fmm_.get_M();
  const auto& T = T_;
  assert(M.M() == T.M());
  assert(M.m() == T.m());
  for (Int br = 0; br < T.M(); ++br) {
    const auto Mb = M.block(br, br);
    const auto Tb = T.block(br, br);
    for (Int k = 0; k < square(M.m()); ++k) {
      diag_num += square(Tb[k] - Mb[k]);
      diag_den += square(Mb[k]);
    }
  }
  printf("> rd(M,T) %1.3e\n", std::sqrt(diag_num/diag_den));
}
