#ifndef INCLUDE_SLMM_MESH_HPP
#define INCLUDE_SLMM_MESH_HPP

#include "slmm_defs.hpp"

namespace slmm {
namespace mesh {

// c is cell (aka element). n is node. e is edge. Hence c2e is the cell-to-edge
// adjacency array.
//   geo is the basic geometric mesh. cgll is a continuous GLL mesh induced by
// the geometric mesh and the reference map. dgll is a discontinuous GLL map
// induced by the CGLL mesh.
//   In a geo cell, the four nodes are ordered CCW. When a geometric mesh is
// converted to GLL, geometric cell i is divided into n = (np-1)^2 subcells. A
// GLL cell is 1-1 with a geo cell and contains GLL subcells.
//   For netcdf I/O, we'll need to make GLL subcells explicitly, and they will
// be numbered n*i : n*(i+1)-1. make_io_cgll_from_internal_cgll does this. We'll
// use 'io' vs 'internal' decoration to distinguish these when necessary.
//   For internal use, we don't need to form these cells explicitly. Instead,
// cgll_c2n has (np-1)^2 slots per slice. Nodes are ordered, e.g. with np=4,
//     12 13 14 15
//      8  9 10 11
//      4  5  6  7
//      0  1  2  3.
// Hence cgll_c2n(i_cell, k) gives the k'th node of cell i_cell.
//   With respect to the reference square (e.g., in siqk::sqr), in a quad the
// bottom-left node is (-1,-1), the bottom-right is (1,0), the top-right is
// (1,1), and the top-left is (-1,1).
//   DGLL topology looks the same except that edge nodes are not
// shared. dglln2cglln(k) maps the k'th DGLL node to the corresponding CGLL
// node. cglln2dglln(k,:) is the list of DGLL nodes associated with CGLL node k.
//   In all topology arrays, -1 indicates the end of a list. E.g., if CGLL node
// k corresponds to 2 DGLL nodes, then cglln2dglln(k,{0,1}) have values, and the
// rest are -1.

void make_cubedsphere(
  Vec3s::HostMirror& geo_p, Idxs::HostMirror& geo_c2n, const Int ne);

void make_cgll_from_geo(
  const Vec3s::HostMirror& geo_p, const Idxs::HostMirror& geo_c2n,
  const Int np, Vec3s::HostMirror& cgll_p, Idxs::HostMirror& cgll_c2n);

void make_io_cgll_from_internal_cgll(
  const Vec3s::HostMirror& cgll_p, const Idxs::HostMirror& cgll_c2n,
  Idxs::HostMirror& cgll_io_c2n);

// dgll_c2n(cell_nmbr, :) contains node numbers for the DGLL mesh. However, a
// separate coordinate array (with redundant coordinates) is not
// created. Instead, use dglln2cglln as follows. The coordinates of the node
// dgll_c2n(cell_nmbr, k) are cgll_p(dglln2cglln(dgll_c2n(cell_nmbr, k)), :).
void make_dgll_from_cgll(
  const Vec3s::HostMirror& cgll_p, const Idxs::HostMirror& cgll_c2n,
  IdxArray::HostMirror& dglln2cglln, Idxs::HostMirror& dgll_c2n);

namespace impl {
// slice(e2n,i) has np slots, and slots 0 and np-1 are filled.
void make_c2e_from_c2n(const Int np, const Idxs::HostMirror& c2n,
                       Idxs::HostMirror& c2e, Idxs::HostMirror& e2n);

// Return 0 if all elements' subtri normals point outward relative to the
// sphere.
Int check_elem_normal_against_sphere(
  const Vec3s::HostMirror& p, const Idxs::HostMirror& e);
} // namespace impl
} // namespace mesh
} // namespace slmm

#endif
