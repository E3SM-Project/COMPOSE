#ifndef INCLUDE_SLMMIR_UTIL_HPP
#define INCLUDE_SLMMIR_UTIL_HPP

#include "slmm_basis.hpp"

#include "slmmir_mesh.hpp"
#include "slmmir.hpp"

void copy_vertices(
  const AVec3s& p, const AIdxs& c2n,
  const Int ci, Real* ps);

Real calc_jacobian(
  const AVec3s& p, const Int* cell, const Real& a,
  const Real& b);

void calc_node_jacobians(
  const Int& np, const Basis& b, const AIdxs& c2n,
  const AVec3s& p, ARealArray& J_dg);

void calc_node_jacobians(
  const Mesh& m, const AVec3s& p, ARealArray& J_dg);

void calc_basis_function_integrals(
  const Int& np, const Int tq_order, const AVec3s& p,
  const AIdxs& c2n, ARealArray& dgbfi);

// Integrals of the basis functions on the sphere(QOS).
void calc_basis_function_integrals(
  const Mesh& m, const AVec3s& p, ARealArray& dgbfi,
  ARealArray& cgbfi);

// Basis function integrals consistent with how Homme defines mass (QOF).
void calc_gll_basis_function_integrals(
  const Int& np, const AIdxs& c2n,
  const AVec3s& p, ARealArray& dgbfi_gll);

void calc_gll_basis_function_integrals(
  const Mesh& m, ARealArray& dgbfi_gll);

void calc_gll_basis_function_integrals(
  const Mesh& m, const Basis& b, ARealArray& dgbfi_gll);

void calc_gll_basis_function_integrals(
  const Mesh& m, const Basis& b, ARealArray& dgbfi, ARealArray& cgbfi);

// The above versions always use the GLL basis weights. This one uses b's.
void calc_gll_basis_function_integrals(
  const Int& np, const Basis& b, const AIdxs& c2n,
  const AVec3s& p, ARealArray& dgbfi_gll);

void calc_gll_basis_function_integrals(
  const Int& np, const Basis& b, const AIdxs& c2n,
  const AVec3s& p, const Real* jacobian, Real* dgbfi_gll);

#endif
