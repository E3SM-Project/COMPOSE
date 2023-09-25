#include "slmm_mesh.hpp"

#include "slmmir_remapper.hpp"
#include "slmmir_util.hpp"

static constexpr Int max_nvert = 8;
static constexpr Int max_hits = 25; // Covers at least a 2-halo.

class CountIntersectionsFunctor {
protected:
  const siqk::sh::Mesh<ko::HostSpace>& cm_;
  const AVec3s p_;
  const AIdxs e_;
  Int hits_[max_hits];
  Int k_, nh_;

public:
  CountIntersectionsFunctor (
    const siqk::sh::Mesh<ko::HostSpace>& cm, const AVec3s& p,
    const AIdxs& c2n)
    : cm_(cm), p_(p), e_(c2n), nh_(0)
  {}

  void reset (const Int clipped_ci) {
    k_ = clipped_ci;
    nh_ = 0;
  }

  void operator() (const Int clip_ci) {
    // Check whether we've clipped against this polygon before and there was a
    // non-0 intersection.
    for (Int i = 0; i < nh_; ++i)
      if (hits_[i] == clip_ci)
        return;
    // We have not, so do the intersection.
    Int no = 0;
    {
      // Area of all overlapping regions.
      // In and out vertex lists.
      Real buf[9*max_nvert];
      siqk::RawVec3s
        vi(buf, max_nvert),
        vo(buf + 3*max_nvert, max_nvert),
        wrk(buf + 6*max_nvert, max_nvert);
      Int ni;
      ni = 0;
      for (Int i = 0; i < szslice(e_); ++i) {
        if (e_(k_,i) == -1) break;
        SphereGeo::copy(slice(vi, i), slice(p_, e_(k_,i)));
        ++ni;
      }
      siqk::sh::clip_against_poly<SphereGeo>(cm_, clip_ci, vi, ni, vo, no, wrk);
    }
    if (no) {
      // Non-0 intersection, so record.
      if (nh_ == max_hits) Kokkos::abort("max_hits is too small.");
      hits_[nh_++] = clip_ci;
    }
  }

  Int get_nhits () const { return nh_; }
  const Int* get_hits () const { return hits_; }
};

static void calc_T_pattern_fwd (
  const Mesh& m, const AVec3s& advected_p,
  const RemapData::Octree& ot,
  Array<Int>& colptr, Array<Int>& rowidx, const bool want_Tgt,
  Array<Int>& rowptr, Array<Int>& colidx)
{
  Timer::start(Timer::ts_remap_T_geometry);
  const Int ncell = nslices(m.geo_c2n);
  AIdxs hits(ncell, max_hits+1);
  {
    siqk::sh::Mesh<ko::HostSpace> cm;
    cm.p = a2ConstVec3s(m.geo_p); cm.e = a2ConstIdxs(m.geo_c2n);
    cm.nml = a2ConstVec3s(m.geo_nml); cm.en = a2ConstIdxs(m.geo_c2nml);
#   pragma omp parallel for schedule(static, 20)
    for (Int ci = 0; ci < ncell; ++ci) {
      Real bb[6];
      RemapData::Octree::calc_bb(advected_p, slice(m.geo_c2n, ci),
                                 szslice(m.geo_c2n), bb);
      CountIntersectionsFunctor cif(cm, advected_p, m.geo_c2n);
      cif.reset(ci);
      ot.apply(bb, cif);
      const Int* ci_hits = cif.get_hits();
      const Int hin = cif.get_nhits();
      hits(ci, 0) = hin;
      for (Int hi = 0; hi < hin; ++hi)
        hits(ci, hi+1) = ci_hits[hi];
    }
  }
  if (want_Tgt) {
    // Make the T' graph, which naturally follows from the above
    // procedure. rowidx is not sorted because we don't need it to be.
    colptr.optclear_and_resize(ncell + 1, 0);
    colptr[0] = 0;
    for (Int ci = 0; ci < ncell; ++ci) {
      const auto hin = hits(ci, 0);
      colptr[ci+1] = colptr[ci] + hin;
    }
    rowidx.optclear_and_resize(colptr[ncell]);
#   pragma omp parallel for
    for (Int ci = 0; ci < ncell; ++ci) {
      Size k = colptr[ci];
      for (Int hi = 1; hi <= hits(ci, 0); ++hi)
        rowidx[k++] = hits(ci, hi);
    }
  }
  Timer::stop(Timer::ts_remap_T_geometry);
  Timer::start(Timer::ts_remap_T_crs);
  // T matrix.
  rowptr.optclear_and_resize(ncell + 1, 0);
  for (Int ci = 0; ci < ncell; ++ci)
    for (Int hi = 1; hi <= hits(ci, 0); ++hi)
      ++rowptr[hits(ci, hi) + 1];
  // Cumsum.
  for (Int ci = 1; ci <= ncell; ++ci)
    rowptr[ci] += rowptr[ci-1];
  colidx.optclear_and_resize(rowptr[ncell]);
  // Shift up 1.
  for (Int ci = ncell; ci > 0; --ci)
    rowptr[ci] = rowptr[ci-1];
  for (Int ci = 0; ci < ncell; ++ci)
    for (Int hi = 1; hi <= hits(ci, 0); ++hi) {
      const Int row = hits(ci, hi);
      colidx[rowptr[row+1]] = ci;
      ++rowptr[row+1];
    }
# pragma omp parallel for
  for (Int ci = 0; ci < ncell; ++ci)
    std::sort(colidx.data() + rowptr[ci], colidx.data() + rowptr[ci+1]);
  Timer::stop(Timer::ts_remap_T_crs);
#ifdef EXPENSIVE_CHECKS
  if (want_Tgt) {
    // Check the transpose.
#   pragma omp parallel for
    for (Int r = 0; r < ncell; ++r)
      for (Int j = rowptr[r]; j < rowptr[r+1]; ++j) {
        const Int c = colidx[j];
        bool found = false;
        for (Int k = colptr[c]; k < colptr[c+1]; ++k)
          if (rowidx[k] == r) {
            found = true;
            break;
          }
        assert(found);
      }
  }
#endif
}

static void partition_n_uniformly (const Int n, const Int nparts,
                                   Array<Int>& p) {
  p.optclear_and_resize(nparts + 1);
  const Int base = n / nparts;
  Int rem = n - base*nparts;
  Int extra = rem > 0 ? 1 : 0;
  p[0] = 0;
  for (Int i = 1; i <= nparts; ++i) {
    p[i] = p[i-1] + base + extra;
    if (rem > 0) {
      --rem;
      if (rem == 0) extra = 0;
    }
  }
}

// Bring in p_s_ol's denominator.
static void finish_p_s_ol (const RemapData::MT& T, Real* const p_s_ol) {
  const Int ncell = T.M();
  const Int np2 = T.m();
  const Size* const rowptr = T.rowptr();
  const Int* const colidx = T.colidx();
  Array<Int> parts;
  partition_n_uniformly(ncell, omp_get_max_threads(), parts);
  Array<Real> p_s_ol_colsum(np2 * ncell);
# pragma omp parallel
  {
    const int tid = omp_get_thread_num();
    const int n = np2*ncell;
#   pragma omp for
    for (Int i = 0; i < n; ++i)
      p_s_ol_colsum[i] = 0;
    for (Int tci = 0; tci < ncell; ++tci) {
      for (Int cj = rowptr[tci]; cj < rowptr[tci+1]; ++cj) {
        const Int sci = colidx[cj];
        if (sci < parts[tid] || sci >= parts[tid+1]) continue;
        Real* p_s_ol_block = p_s_ol + np2*cj;
        Real* const p_s_ol_colsum_block = p_s_ol_colsum.data() + sci*np2;
        for (Int i = 0; i < np2; ++i)
          p_s_ol_colsum_block[i] += p_s_ol_block[i];
      }
    }
#   pragma omp barrier
#   pragma omp for schedule(static, 1)
    for (Int tci = 0; tci < ncell; ++tci) {
      Real* p_s_ol_block = p_s_ol + np2 * rowptr[tci];
      for (Int cj = rowptr[tci]; cj < rowptr[tci+1]; ++cj) {
        const Int sci = colidx[cj];
        Real* const p_s_ol_colsum_block = p_s_ol_colsum.data() + sci*np2;
        for (Int i = 0; i < np2; ++i)
          p_s_ol_block[i] /= p_s_ol_colsum_block[i];
        p_s_ol_block += np2;
      }
    }
  }
}

static void sum_cols (const RemapData::MT& T, Real* const colsum) {
  const Int ncell = T.M();
  const Int np2 = T.m(), np4 = np2*np2;
  const Size* const rowptr = T.rowptr();
  const Int* const colidx = T.colidx();
  Array<Int> parts;
  partition_n_uniformly(ncell, omp_get_max_threads(), parts);
# pragma omp parallel
  {
    const int tid = omp_get_thread_num();
    const int n = np2*ncell;
#   pragma omp for
    for (Int i = 0; i < n; ++i)
      colsum[i] = 0;
    for (Int tci = 0; tci < ncell; ++tci) {
      const Int cj0 = rowptr[tci];
      for (Int cj = cj0; cj < rowptr[tci+1]; ++cj) {
        const Int sci = colidx[cj];
        if (sci < parts[tid] || sci >= parts[tid+1]) continue;
        const Real* block = T.blockrow(tci) + np4*(cj - cj0);
        Real* const colsum_block = colsum + sci*np2;
        for (Int j = 0; j < np2; ++j) {
          // Little block column sum.
          Real sum = 0;
          for (Int i = 0; i < np2; ++i)
            sum += block[i*np2 + j];
          // Accum into big block column sum.
          colsum_block[j] += sum;
        }
      }
    }
  }
}

static void fill_T_fwd_facet (
  const Mesh& m, const AVec3s& advected_p, RemapData::MT& T,
  const Method::Enum method, Array<Real>* p_s_ol)
{
  const Int ncell = nslices(m.geo_c2n);
  const Int np = m.np, np2 = square(np), np4 = square(np2);
  siqk::TriangleQuadrature tq;
  siqk::RawConstVec3s tq_bary;
  siqk::RawConstArray tq_w;
  tq.get_coef(m.tq_order, tq_bary, tq_w);
  const Int nq = len(tq_w);
  GLL gll;
  const Size* rowptr = T.rowptr();
  const Int* colidx = T.colidx();
  siqk::sh::Mesh<ko::HostSpace> cm;
  cm.p = a2ConstVec3s(m.geo_p); cm.e = a2ConstIdxs(m.geo_c2n);
  cm.nml = a2ConstVec3s(m.geo_nml); cm.en = a2ConstIdxs(m.geo_c2nml);  
  if (p_s_ol) p_s_ol->optclear_and_resize(np2 * rowptr[ncell]);
# pragma omp parallel for schedule(static, 1)
  for (Int tci = 0; tci < ncell; ++tci) {
    Real* block = T.blockrow(tci);
    Real* p_s_ol_block = p_s_ol ? p_s_ol->data() + np2 * rowptr[tci] : nullptr;
    const auto tcell = slice(m.geo_c2n, tci);
    for (Int cj = rowptr[tci]; cj < rowptr[tci+1]; ++cj) {
      const Int sci = colidx[cj];
      const auto scell = slice(m.geo_c2n, sci);
      Real buf[9*max_nvert];
      // Get overlap element.
      Int no;
      siqk::RawVec3s vo(buf, max_nvert); {
        siqk::RawVec3s
          vi(buf + 3*max_nvert, max_nvert),
          wrk(buf + 6*max_nvert, max_nvert);
        Int ni = 0;
        for (Int i = 0; i < szslice(m.geo_c2n); ++i) {
          if (scell[i] == -1) break;
          SphereGeo::copy(slice(vi, i), slice(advected_p, scell[i]));
          ++ni;
        }
        siqk::sh::clip_against_poly<SphereGeo>(cm, tci, vi, ni, vo, no, wrk);
        assert(no);
      }
      // Get overlap element vertices in ref elements.
      Real* const tvo = buf + 3*no;
      Real* const svo = tvo + 2*no;
      for (Int k = 0; k < no; ++k) {
        // For reasons that are not clear to me, the global mass is sensitive to
        // the quality of the following nonlinear solves. In one way, that's
        // obvious: one needs a watertight and nonoverlapping covering of the
        // faceted object to get mass conservation to machine precision. But the
        // puzzle is why QOS's similar calculations (in fact, at every
        // quadrature point) don't need more accuracy. Perhaps somehow the
        // faceted shape's area is harder to compute.
        static const Real tol = std::numeric_limits<Real>::epsilon();
        static const int max_nits = 10;
        siqk::sqr::calc_sphere_to_ref(m.geo_p, tcell, slice(vo,k),
                                      tvo[2*k], tvo[2*k+1],
                                      nullptr, max_nits, tol);
        siqk::sqr::calc_sphere_to_ref(advected_p, scell, slice(vo,k),
                                      svo[2*k], svo[2*k+1],
                                      nullptr, max_nits, tol);
      }
      for (Int i = 0; i < np4; ++i) block[i] = 0;
      if (p_s_ol)
        for (Int i = 0; i < np2; ++i) p_s_ol_block[i] = 0;
      for (Int ktri = 1; ktri < no-1; ++ktri) { // triangles in vo
        const Real tri_area =
          Method::is_ir(method)
          // Here, integration is over x, with no J_x^a factor, because this is
          // just a projection; no density factor enters this part of the
          // calculation.
          ? PlaneGeo::calc_tri_jacobian(tvo, tvo+2*ktri, tvo+2*(ktri+1))
          // J_x^a dx = J_reftri^a da. Thus, for QOF CDG, integration is done
          // over a, whereas for QOS CDG, it is done over x.
          : PlaneGeo::calc_tri_jacobian(svo, svo+2*ktri, svo+2*(ktri+1));
        const Real f = 0.5 * tri_area;
        for (Int q = 0; q < nq; ++q) { // quad point
          Real tgj[GLL::np_max], tgi[GLL::np_max], sgj[GLL::np_max],
            sgi[GLL::np_max];
          {
            Real tvo_coord[2];
            PlaneGeo::bary2coord(tvo, tvo+2*ktri, tvo+2*(ktri+1),
                                 slice(tq_bary, q), tvo_coord);
            gll.eval(m.np, tvo_coord[0], tgi);
            gll.eval(m.np, tvo_coord[1], tgj);
          }
          {
            Real svo_coord[2];
            PlaneGeo::bary2coord(svo, svo+2*ktri, svo+2*(ktri+1),
                                 slice(tq_bary, q), svo_coord);
            gll.eval(m.np, svo_coord[0], sgi);
            gll.eval(m.np, svo_coord[1], sgj);
          }
          Real d0 = f * tq_w[q];
          if (p_s_ol) {
            for (Int sj = 0, s_basis_idx = 0; sj < np; ++sj) {
              const Real d1 = d0 * sgj[sj];
              for (Int si = 0; si < np; ++si, ++s_basis_idx) {
                const Real d = d1 * sgi[si];
                p_s_ol_block[s_basis_idx] += d;
              }
            }
          }
          for (Int tj = 0, t_basis_idx = 0; tj < np; ++tj) {
            const Real d1 = d0 * tgj[tj];
            for (Int ti = 0; ti < np; ++ti, ++t_basis_idx) {
              Real d2 = d1 * tgi[ti];
              for (Int sj = 0, s_basis_idx = 0; sj < np; ++sj) {
                const Real d3 = d2 * sgj[sj];
                for (Int si = 0; si < np; ++si, ++s_basis_idx) {
                  const Real d = d3 * sgi[si];
                  block[np2*t_basis_idx + s_basis_idx] += d;
                }
              }
            }
          }
        }
      }
      block += np4;
      if (p_s_ol) p_s_ol_block += np2;
    }
  }
  if (p_s_ol) finish_p_s_ol(T, p_s_ol->data());
}

static void fill_T_fwd_sphere (
  const Mesh& m, const AVec3s& advected_p, RemapData::MT& T,
  const Method::Enum method, Array<Real>* p_s_ol)
{
  const Int ncell = nslices(m.geo_c2n);
  const Int np = m.np, np2 = square(np), np4 = square(np2);
  siqk::TriangleQuadrature tq;
  siqk::RawConstVec3s tq_bary;
  siqk::RawConstArray tq_w;
  tq.get_coef(m.tq_order, tq_bary, tq_w);
  const Int nq = len(tq_w);
  GLL gll;
  const Size* rowptr = T.rowptr();
  const Int* colidx = T.colidx();
  siqk::sh::Mesh<ko::HostSpace> cm;
  cm.p = a2ConstVec3s(m.geo_p); cm.e = a2ConstIdxs(m.geo_c2n);
  cm.nml = a2ConstVec3s(m.geo_nml); cm.en = a2ConstIdxs(m.geo_c2nml);  
  if (p_s_ol) p_s_ol->optclear_and_resize(np2 * rowptr[ncell]);
# pragma omp parallel for schedule(static, 1)
  for (Int tci = 0; tci < ncell; ++tci) {
    Real* block = T.blockrow(tci);
    Real* p_s_ol_block = p_s_ol ? p_s_ol->data() + np2 * rowptr[tci] : nullptr;
    const auto tcell = slice(m.geo_c2n, tci);
    for (Int cj = rowptr[tci]; cj < rowptr[tci+1]; ++cj) {
      const Int sci = colidx[cj];
      const auto scell = slice(m.geo_c2n, sci);
      Real buf[9*max_nvert];
      siqk::RawVec3s
        vi(buf, max_nvert),
        vo(buf + 3*max_nvert, max_nvert),
        wrk(buf + 6*max_nvert, max_nvert);
      Int ni = 0, no;
      for (Int i = 0; i < szslice(m.geo_c2n); ++i) {
        if (scell[i] == -1) break;
        SphereGeo::copy(slice(vi, i), slice(advected_p, scell[i]));
        ++ni;
      }
      siqk::sh::clip_against_poly<SphereGeo>(cm, tci, vi, ni, vo, no, wrk);
      assert(no);
      for (Int i = 0; i < np4; ++i) block[i] = 0;
      if (p_s_ol)
        for (Int i = 0; i < np2; ++i) p_s_ol_block[i] = 0;
      for (Int ktri = 1; ktri < no-1; ++ktri) { // triangles in vo
        for (Int q = 0; q < nq; ++q) { // quad point
          Real sphere_coord[3];
          const Real jac_reftri2sphere = SphereGeo::calc_tri_jacobian(
            slice(vo,0), slice(vo,ktri), slice(vo,ktri+1), slice(tq_bary,q),
            sphere_coord);
          Real d0 = 0.5 * tq_w[q] * jac_reftri2sphere;
          Real tgj[GLL::np_max], tgi[GLL::np_max], sgj[GLL::np_max],
            sgi[GLL::np_max];
          {
            Real ta, tb, sa, sb;
            siqk::sqr::calc_sphere_to_ref(m.geo_p   , tcell, sphere_coord,
                                          ta, tb);
            siqk::sqr::calc_sphere_to_ref(advected_p, scell, sphere_coord,
                                          sa, sb);
            gll.eval(m.np, tb, tgj);
            gll.eval(m.np, ta, tgi);
            gll.eval(m.np, sb, sgj);
            gll.eval(m.np, sa, sgi);
            if ( ! Method::is_ir(method))
              d0 *= (calc_jacobian(m.geo_p   , scell, sa, sb) /
                     calc_jacobian(advected_p, scell, sa, sb));
          }
          if (p_s_ol) {
            for (Int sj = 0, s_basis_idx = 0; sj < np; ++sj) {
              const Real d1 = d0 * sgj[sj];
              for (Int si = 0; si < np; ++si, ++s_basis_idx) {
                const Real d = d1 * sgi[si];
                p_s_ol_block[s_basis_idx] += d;
              }
            }
          }
          for (Int tj = 0, t_basis_idx = 0; tj < np; ++tj) {
            const Real d1 = d0 * tgj[tj];
            for (Int ti = 0; ti < np; ++ti, ++t_basis_idx) {
              Real d2 = d1 * tgi[ti];
              for (Int sj = 0, s_basis_idx = 0; sj < np; ++sj) {
                const Real d3 = d2 * sgj[sj];
                for (Int si = 0; si < np; ++si, ++s_basis_idx) {
                  const Real d = d3 * sgi[si];
                  block[np2*t_basis_idx + s_basis_idx] += d;
                }
              }
            }
          }
        }
      }
      block += np4;
      if (p_s_ol) p_s_ol_block += np2;
    }
  }
  if (p_s_ol) finish_p_s_ol(T, p_s_ol->data());
}

static void calc_T_fwd (const Mesh& m, const AVec3s& advected_p,
                        RemapData& rd, const bool want_Tgt)
{
  { // Build T's sparse matrix nonzero pattern.
    Array<Int> rowptr, colidx;
    calc_T_pattern_fwd(m, advected_p, rd.octree(), rd.Tgt().colptr,
                       rd.Tgt().rowidx, want_Tgt, rowptr, colidx);
    const Int N = len(rowptr)-1, n = square(m.np);
    rd.T().init(N, N, n, n, rowptr.data(), colidx.data());
  }
  Timer::start(Timer::ts_remap_T_fill);
  const auto p_s_ol = (Dmc::is_locally_constrained(rd.dmc()) ?
                       &rd.p_s_ol() : nullptr);
  if (Dmc::is_facet(rd.dmc()))
    fill_T_fwd_facet(m, advected_p, rd.T(), rd.method(), p_s_ol);
  else
    fill_T_fwd_sphere(m, advected_p, rd.T(), rd.method(), p_s_ol);
  Timer::stop(Timer::ts_remap_T_fill);
}

void Remapper::
dgll_set_src_tgt_in (Real* const src_in, Real* const tgt_in,
                     Real*& src, Real*& tgt) {
  src = src_in;
  tgt = tgt_in;
}

void Remapper::
dgll_set_tgt_out (const Real* const tgt, Real* const tgt_out) {
}

void Remapper::
project_nolimiter (
  Real* const src_tracer, Real* const tgt_tracer, const Int ntracers,
  Real* const src_density, Real* const tgt_density)
{
  const Int len = dnn_;
  {
    Real* src, * tgt;
    dgll_set_src_tgt_in(src_density, tgt_density, src, tgt);
    apply_R(src, tgt);
    dgll_set_tgt_out(tgt, tgt_density);
  }
  for (Int ti = 0; ti < ntracers; ++ti) {
    Real* src, * tgt;
    dgll_set_src_tgt_in(src_tracer + ti*len, tgt_tracer + ti*len, src, tgt);
    apply_R(src, tgt);
    dgll_set_tgt_out(tgt, tgt_tracer + ti*len);
  }
}

void Remapper::
apply_R (const Real* const src, Real* const tgt) const {
  rd_->remap(src, dnn_, tgt, dnn_, 1, FsmoFtm_.data());
}

void Remapper::
apply_R (const Real* const src, Real* const tgt,
         const Int* const tis, const Int nti) const {
  if ( ! nti) {
    apply_R(src, tgt);
    return;
  }
# pragma omp parallel for
  for (Int i = 0; i < nti; ++i) {
    const Int ti = tis[i];
    rd_->remap_cell(ti, src, dnn_, tgt, dnn_, 1, FsmoFtm_.data());
  }
}

void Remapper::
perturb_rho (Real* const rho, const Real p) {
  const auto& F = rd_->dgbfi_mass();
  const Int n = F.size();
  static Array<Real> u(n);
  Real rho_min = rho[0], utF = 0, FtF = 0, u_max = 0;
# pragma omp parallel
  {
#   pragma omp for reduction (min: rho_min)
    for (Int i = 0; i < n; ++i) rho_min = std::min(rho_min, rho[i]);
    // Random numbers.
#   pragma omp for
    for (Int i = 0; i < n; ++i) u[i] = 2*urand() - 1;
    // Project out the global mass.
#   pragma omp for reduction (+: utF)
    for (Int i = 0; i < n; ++i) utF += F[i]*u[i];
#   pragma omp for reduction (+: FtF)
    for (Int i = 0; i < n; ++i) FtF += F[i]*F[i];
    Real f = utF/FtF;
#   pragma omp for
    for (Int i = 0; i < n; ++i) u[i] -= f*F[i];
    // Scale so that perturbed rho is >= 0.
#   pragma omp for reduction (max: u_max)
    for (Int i = 0; i < n; ++i) u_max = std::max(u_max, std::abs(u[i]));
    f = p*rho_min/u_max;
    // Perturb rho.
#   pragma omp for
    for (Int i = 0; i < n; ++i) rho[i] += f*u[i];
  }
}

void Remapper::
project_and_limit_cdr (
  Real* const src_tracer, Real* const tgt_tracer, const Int ntracers,
  Real* const src_density, Real* const tgt_density,
  const Filter::Enum cdr_method, const Real rho_perturbation)
{
  if ( ! mrd_) {
    const Int ne = std::sqrt(ncell_ / 6);
    mrd_ = std::make_shared<spf::MassRedistributor>(
      ne, rd_->dgbfi_mass(), Filter::to_mrd(cdr_method));
  }
    
  const Int len = dnn_;

  if (record_total_mass_redistribution_)
    total_mass_redistribution_.resize(ntracers+1);

  { // Density.
    Real* src, * tgt;
    dgll_set_src_tgt_in(src_density, tgt_density, src, tgt);
    apply_R(src, tgt);
    // Iterate over target cells and limit the target field.
#   pragma omp parallel for
    for (Int ci = 0; ci < ncell_; ++ci) {
      //const Int info =
      md_->limit_density(ci, tgt + ci*np2_);
      // Technically, this QP can fail. It can be fixed by projecting the
      // source rho to a monotone basis. But in practice this is not an issue
      // for atmospheric transport, since density doesn't get near 0 relative
      // to maximum density.
      //assert(info >= 0);
    }
    dgll_set_tgt_out(tgt, tgt_density);
    if (rho_perturbation)
      perturb_rho(tgt_density, rho_perturbation);
    copy(rho_src_cmbc_.data(), src, dnn_);
    copy(rho_tgt_cmbc_.data(), tgt, dnn_);
  }

  // Tracers.
  for (Int tri = 0; tri < ntracers; ++tri) {
    Real* src, * tgt;
    dgll_set_src_tgt_in(src_tracer + tri*len, tgt_tracer + tri*len, src, tgt);
    // Remap and record cell mass bounds and discrepancy.
#   pragma omp parallel for
    for (Int ti = 0; ti < ncell_; ++ti) {
      rd_->remap_cell(ti, src, dnn_, tgt, dnn_, 1, FsmoFtm_.data());
      Real q_min, q_max;
      md_->calc_q_min_max(rd_->T(), ti, rho_src_cmbc_.data(), src,
                          q_min, q_max);
      mrd_->record(ti, q_min, q_max, rho_tgt_cmbc_.data(), tgt);
    }
    // Quasi-local mass redistribution based on cell mean boundedness.
    mrd_->redistribute_mass();
    if (record_total_mass_redistribution_) {
      Real tmr = 0;
#     pragma omp parallel for reduction(+: tmr)
      for (Int i = 0; i < ncell_; ++i)
        tmr += std::abs(mrd_->get_delta_mass(i));
      total_mass_redistribution_[tri] = tmr;
    }
    // Shape preservation filter with cell mass deltas.
#   pragma omp parallel for
    for (Int ti = 0; ti < ncell_; ++ti) {
      const Real q_min = mrd_->get_q_min(ti);
      const Real q_max = mrd_->get_q_max(ti);
      const Real Q_extra = mrd_->get_delta_mass(ti);
      Int info = md_->limit_tracer(
        ti, q_min, q_max, rho_tgt_cmbc_.data() + ti*np2_, tgt + ti*np2_,
        false, Q_extra);
      if (info < 0) {
        // FP means the cell mean isn't exactly bounded. Solve the QP to get
        // close.
        info = md_->limit_tracer(
          ti, q_min, q_max, rho_tgt_cmbc_.data() + ti*np2_, tgt + ti*np2_,
          true, Q_extra);
        assert_expensive(info >= 0);
      }
    }
    dgll_set_tgt_out(tgt, tgt_tracer + tri*len);
  }    
}

void Remapper::
project_and_limit (
  Real* const src_tracer, Real* const tgt_tracer, const Int ntracers,
  Real* const src_density, Real* const tgt_density)
{
  switch (md_->global_type()) {
  case Filter::qlt:
  case Filter::caas:
  case Filter::mn2:
    project_and_limit_cdr(src_tracer, tgt_tracer, ntracers,
                          src_density, tgt_density,
                          md_->global_type(), perturb_rho_);
    break;
  default:
    throw std::runtime_error("md_->type() is not valid.");
  }
}

Remapper
::Remapper (const PRefineData::Ptr& pr, const MonoData::Ptr& md, const D2Cer::Ptr& d2cer)
  : np_(pr->m_t->np), np2_(np_*np_), ncell_(nslices(pr->m_t->geo_c2n)),
    dnn_(ncell_*np2_), cnn_(nslices(pr->m_t->cgll_p)), m_(pr->m_t), rd_(pr->rd_t),
    md_(md), d2cer_(d2cer), perturb_rho_(0), record_total_mass_redistribution_(false),
    subcell_bounds_(false), pr_(pr)
{
  if (md_) {
    resize(rho_src_cmbc_, dnn_);
    resize(rho_tgt_cmbc_, dnn_);
  }
  if (Method::is_isl(rd_->method()))
    init_isl();
}

void Remapper::
use_subcell_bounds () {
  SIQK_THROW_IF( ! Method::is_isl(rd_->method()),
                 "Subcell bounds are for ISL only.");
  subcell_bounds_ = true;
  resize(q_data_, square(np_-1)*ncell_, 2);
}

void Remapper::
use_fit_extremum () {
  SIQK_THROW_IF( ! Method::is_isl(rd_->method()),
                 "Quadratic extremum is impl'ed ISL only.");
  fit_extremum_ = std::make_shared<FitExtremum>(np_);
}

// Stuff to analyze tracer consistency via project_and_limit_cdr.
void Remapper::
set_rho_perturbation (const Real p) { perturb_rho_ = p; }

void Remapper::
record_total_mass_redistribution (const bool record) {
  record_total_mass_redistribution_ = record;
}  

Real Remapper::
get_total_redistributed_mass (const Int& tracer_idx) {
  if ( ! record_total_mass_redistribution_) return 0;
  return total_mass_redistribution_[tracer_idx];
}

Real Remapper::
get_total_mass_discrepancy (const Int& tracer_idx) {
  if ( ! record_total_mass_redistribution_) return 0;
  return total_mass_discrepancy_.empty() ? 0 :
    total_mass_discrepancy_[tracer_idx];
}

// On input, src_tracer is rho*tracer. On output, it is just the updated
// tracer. Density is removed for output and error checking.
void Remapper::
remap (const AVec3s& advected_p,
       Real* const src_tracer, Real* const tgt_tracer, const Int ntracers,
       Real* const src_density, Real* const tgt_density)
{
  // Compute T and p_s_ol.
  Timer::start(Timer::ts_remap_T);
  calc_T_fwd(*m_, advected_p, *rd_, md_ != nullptr);
  Timer::stop(Timer::ts_remap_T);

  if (Method::is_ir(rd_->method())) {
    // Compute density factor.
    Timer::start(Timer::ts_remap_node_jac);
    if (Dmc::is_facet(rd_->dmc())) {
      resize(FsmoFtm_, rd_->Jt().size());
      // Advected source basis function eta_ij's integral is the sum of its
      // mixed mass matrix column.
      sum_cols(rd_->T(), FsmoFtm_.data());
      // Eulerian source basis function eta_ij's integral is the GLL weight.
      GLL gll;
      const Real* gll_x, * gll_wt;
      gll.get_coef(m_->np, gll_x, gll_wt);
      const Int np = m_->np;
#     pragma omp parallel for
      for (Int si = 0; si < ncell_; ++si) {
        Real* const FsmoFtm = &FsmoFtm_[np2_*si];
        for (Int i = 0, k = 0; i < np; ++i)
          for (Int j = 0; j < np; ++j, ++k)
            FsmoFtm[k] = (gll_wt[i]*gll_wt[j])/FsmoFtm[k];
      }
    } else {
      calc_basis_function_integrals(m_->np, m_->tq_order, advected_p,
                                    m_->geo_c2n, FsmoFtm_);
#     pragma omp parallel for
      for (Int i = 0; i < dnn_; ++i)
        FsmoFtm_[i] = rd_->dgbfi()[i]/FsmoFtm_[i];
    }
    Timer::stop(Timer::ts_remap_node_jac);
  }

  Timer::start(Timer::ts_remap_project);
  if ( ! md_)
    project_nolimiter(src_tracer, tgt_tracer, ntracers,
                      src_density, tgt_density);
  else
    project_and_limit(src_tracer, tgt_tracer, ntracers,
                      src_density, tgt_density);
  Timer::stop(Timer::ts_remap_project);
}
