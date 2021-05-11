#include "slmmir_util.hpp"
#include "slmm_gll.hpp"

void copy_vertices (
  const AVec3s& p, const AIdxs& c2n,
  const Int ci, Real* ps)
{
  const auto cell = slice(c2n, ci);
  for (Int i = 0; i < szslice(c2n); ++i) {
    const auto n = slice(p, cell[i]);
    for (Int k = 0; k < 3; ++k) ps[k] = n[k];
    ps += 3;
  }
}

Real calc_jacobian (
  const AVec3s& p, const Int* cell, const Real& a,
  const Real& b)
{
  Real J[9];
  siqk::sqr::impl::calc_Jacobian(p, cell, a, b, J);
  SphereGeo::cross(J, J+3, J+6);
  return std::sqrt(SphereGeo::norm2(J+6));
}

void calc_node_jacobians (
  const Int& np, const Basis& b, const AIdxs& c2n,
  const AVec3s& p, ARealArray& J_dg)
{
  const Int np2 = square(np);
  resize(J_dg, nslices(c2n)*np2);
  const Real* xnode;
  b.get_x(np, xnode);
# pragma omp parallel for
  for (Int ci = 0; ci < nslices(c2n); ++ci) {
    const auto cell = slice(c2n, ci);
    for (Int j = 0, basis_idx = 0; j < np; ++j)
      for (Int i = 0; i < np; ++i, ++basis_idx)
        J_dg(ci*np2 + basis_idx) = calc_jacobian(p, cell, xnode[i], xnode[j]);
  }  
}

void calc_node_jacobians (
  const Mesh& m, const AVec3s& p, ARealArray& J_dg)
{
  calc_node_jacobians(m.np, *m.basis, m.geo_c2n, p, J_dg);
}

void calc_basis_function_integrals (
  const Int& np, const Int tq_order, const AVec3s& p,
  const AIdxs& c2n, ARealArray& dgbfi)
{
  const Int np2 = square(np);
  resize(dgbfi, nslices(c2n)*np2);
  copy(dgbfi, 0);
  siqk::TriangleQuadrature tq;
  siqk::RawConstVec3s tq_bary;
  siqk::RawConstArray tq_w;
  tq.get_coef(tq_order, tq_bary, tq_w);
  const Int nq = len(tq_w);
  GLL b;
# pragma omp parallel for
  for (Int ci = 0; ci < nslices(c2n); ++ci) { // cell
    Real ps[12];
    copy_vertices(p, c2n, ci, ps);
    const auto cell = slice(c2n, ci);
    for (Int k = 1; k <= 2; ++k) // 2 triangles per quad cell
      for (Int q = 0; q < nq; ++q) { // quad point
        Real sphere_coord[3];
        const Real jac = SphereGeo::calc_tri_jacobian(
          ps, ps+3*k, ps+3*(k+1), slice(tq_bary, q), sphere_coord);
        Real gj[GLL::np_max], gi[GLL::np_max]; {
          Real al, be;
          siqk::sqr::calc_sphere_to_ref(p, cell, sphere_coord, al, be);
          b.eval(np, be, gj);
          b.eval(np, al, gi);
        }
        const Real d0 = 0.5 * tq_w[q] * jac;
        for (Int j = 0, basis_idx = 0; j < np; ++j) { // along ref y dir
          const Real d1 = d0 * gj[j];
          for (Int i = 0; i < np; ++i, ++basis_idx) // along ref x dir
            dgbfi(ci*np2 + basis_idx) += d1 * gi[i];
        }
      }
  }
}

// Integrals of the basis functions on the sphere (QOS).
void calc_basis_function_integrals (
  const Mesh& m, const AVec3s& p, ARealArray& dgbfi,
  ARealArray& cgbfi)
{
  calc_basis_function_integrals(m.np, m.tq_order, p, m.geo_c2n, dgbfi);
  resize(cgbfi, nslices(m.cgll_p));
  copy(cgbfi, 0);
  for (Int i = 0; i < len(m.dglln2cglln); ++i)
    cgbfi(m.dglln2cglln(i)) += dgbfi(i);
}

// Basis function integrals consistent with how Homme defines mass (QOF).
void calc_gll_basis_function_integrals (
  const Int& np, const Basis& b, const AIdxs& c2n,
  const AVec3s& p, ARealArray& dgbfi_gll)
{
  const Int np2 = square(np);
  resize(dgbfi_gll, nslices(c2n)*np2);
  const Real* xnode, * wt;
  b.get_x(np, xnode);
  b.get_w(np, wt);
  {
    bool fnd = false;
    for (int i = 0; i < np; ++i) if (wt[i] <= 0) fnd = true;
    if (fnd) std::cerr << "WARNING: some w < 0\n";
  }
# pragma omp parallel for
  for (Int ci = 0; ci < nslices(c2n); ++ci) {
    const auto cell = slice(c2n, ci);
    for (Int j = 0, basis_idx = 0; j < np; ++j)
      for (Int i = 0; i < np; ++i, ++basis_idx) {
        const Real jac = calc_jacobian(p, cell, xnode[i], xnode[j]);
        // Product of weights is the integral of the 2D basis function on the
        // ref square. Multiply by Jacobian of the map bilinear quad ->
        // sphere. Since this is GLL quadrature, there's exactly one quadrature
        // point.
        dgbfi_gll(ci*np2 + basis_idx) = 0.25 * jac * wt[i] * wt[j];
      }
  }
}

void calc_gll_basis_function_integrals (
  const Int& np, const AIdxs& c2n,
  const AVec3s& p, ARealArray& dgbfi_gll)
{
  calc_gll_basis_function_integrals(np, GLL(), c2n, p, dgbfi_gll);
}

void calc_gll_basis_function_integrals (
  const Mesh& m, ARealArray& dgbfi_gll)
{
  calc_gll_basis_function_integrals(m.np, *m.basis, m.geo_c2n, m.geo_p, dgbfi_gll);
}

void calc_gll_basis_function_integrals (
  const Mesh& m, const Basis& b, ARealArray& dgbfi_gll)
{
  calc_gll_basis_function_integrals(m.np, b, m.geo_c2n, m.geo_p, dgbfi_gll);
}

void calc_gll_basis_function_integrals (
  const Mesh& m, const Basis& b, ARealArray& dgbfi, ARealArray& cgbfi)
{
  calc_gll_basis_function_integrals(m, b, dgbfi);
  resize(cgbfi, nslices(m.cgll_p));
  copy(cgbfi, 0);
  for (Int i = 0; i < len(m.dglln2cglln); ++i)
    cgbfi(m.dglln2cglln(i)) += dgbfi(i);
}

void calc_gll_basis_function_integrals (
  const Int& np, const Basis& b, const AIdxs& c2n,
  const AVec3s& p, const Real* jacobian, Real* dgbfi_gll)
{
  const Int np2 = square(np);
  const Real* wt;
  b.get_w(np, wt);
  {
    bool fnd = false;
    for (int i = 0; i < np; ++i) if (wt[i] <= 0) fnd = true;
    if (fnd) std::cerr << "WARNING: some w < 0\n";
  }
# pragma omp parallel for
  for (Int ci = 0; ci < nslices(c2n); ++ci) {
    for (Int j = 0, basis_idx = 0; j < np; ++j)
      for (Int i = 0; i < np; ++i, ++basis_idx) {
        const Int idx = ci*np2 + basis_idx;
        dgbfi_gll[idx] = 0.25 * jacobian[idx] * wt[i] * wt[j];
      }
  }
}
