#ifndef INCLUDE_SLMM_MESH_HPP
#define INCLUDE_SLMM_MESH_HPP

#include "slmm_array_tree.hpp"
#include "slmm_mesh.hpp"

namespace slmm {

struct Basis;

namespace mesh {

// c is cell (aka element). n is node. e is edge. Hence c2e is the cell-to-edge
// adjacency array. sc is subcell.
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
// cgll_c2n has np^2 slots per slice. Nodes are ordered, e.g. with np=4,
//     12 13 14 15
//      8  9 10 11
//      4  5  6  7
//      0  1  2  3.
// Hence cgll_c2n(i_cell, k) gives the k'th node of cell i_cell.
//   With respect to the reference square (e.g., in siqk::sqr), in a quad the
// bottom-left node is (-1,-1), the bottom-right is (1,-1), the top-right is
// (1,1), and the top-left is (-1,1).
//   DGLL topology looks the same except that edge nodes are not
// shared. dglln2cglln(k) maps the k'th DGLL node to the corresponding CGLL
// node.
//   In all topology arrays, -1 indicates the end of a list. E.g., if CGLL node
// k corresponds to 2 DGLL nodes, then cglln2dglln(k,{0,1}) have values, and the
// rest are -1.

void make_cubedsphere_mesh(
  AVec3s& geo_p, AIdxs& geo_c2n, const Int ne);

void make_nonuniform(AVec3s& geo_p);

// Rotate the grid given (axis, angle). Store the rotation in the row-major 3x3
// matrix R.
void rotate_grid(const Real axis[3], const Real angle, Real* R, AVec3s& p);

// Make geometric data for a mesh that is a cubed-sphere refined into a subcell
// mesh, by default the CGLL mesh. This is used to study p=1 transport over a
// p>1 GLL mesh, for example. In that case, call this instead of
// make_cubedsphere_mesh.
void make_cubedsphere_subcell_mesh(
  const Int ne, const Int np, const Basis& basis,
  AVec3s& geo_p, AIdxs& geo_c2n);

void make_cgll_from_geo(
  const AVec3s& geo_p, const AIdxs& geo_c2n,
  const Int np, const Basis& basis, AVec3s& cgll_p,
  AIdxs& cgll_c2n);

// ref_x is the [-1,1] reference coords for the element.
void make_subcell_from_geo(
  const AVec3s& geo_p, const AIdxs& geo_c2n, const Int np,
  const Real* const ref_x, AVec3s& sc_p, AIdxs& sc_c2n);

void make_io_cgll_from_internal_cgll(
  const AVec3s& cgll_p, const AIdxs& cgll_c2n,
  AIdxs& cgll_io_c2n);

// dgll_c2n(cell_nmbr, :) contains node numbers for the DGLL mesh. However, a
// separate coordinate array (with redundant coordinates) is not
// created. Instead, use dglln2cglln as follows. The coordinates of the node
// dgll_c2n(cell_nmbr, k) are cgll_p(dglln2cglln(dgll_c2n(cell_nmbr, k)), :).
void make_dgll_from_cgll(
  const AVec3s& cgll_p, const AIdxs& cgll_c2n,
  AIdxArray& dglln2cglln, AIdxs& dgll_c2n);

// Compute a compressed array map of a cell to its cell neighbors.
void get_adjacent_cells(const AIdxs& geo_c2n, AIdxArray& geo_c2cnbrs_ptr,
                        AIdxArray& geo_c2cnbrs);

// Given a coordinate (x,y,z) on the sphere, return the index of the cell it is
// in for the quasiuniform cubed-sphere mesh with ne x ne faces.
//   If angle != 0, then the mesh was rotated by row-major rotation matrix R.
Int get_cell_idx(const Int ne, const Real angle, const Real* const R,
                 Real x, Real y, Real z);

// tree is a tree over geometric cells. It is stored as an Int array.
bool make_cubedsphere_tree_over_cells(const Int ne, AIdxArray& tree,
                                      const Int n_user_slots=0,
                                      const bool perform_checks=false);

namespace impl {
// slice(e2n,i) has np slots, and slots 0 and np-1 are filled.
void make_c2e_from_c2n(const Int np, const AIdxs& c2n,
                       AIdxs& c2e, AIdxs& e2n);

// Return 0 if all elements' subtri normals point outward relative to the
// sphere.
Int check_elem_normal_against_sphere(
  const AVec3s& p, const AIdxs& e);
} // namespace impl
} // namespace mesh
} // namespace slmm

#endif
