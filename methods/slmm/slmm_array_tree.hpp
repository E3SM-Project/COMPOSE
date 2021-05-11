#ifndef INCLUDE_SLMM_ARRAY_TREE_HPP
#define INCLUDE_SLMM_ARRAY_TREE_HPP

#include "slmm_defs.hpp"

namespace slmm {
namespace tree {

// This is an array representation of a tree. An entry in this array has one of
// five types. They are arranged in this order:
//   * number of nodes or cells in the node;
//   * parent node index;
//   * negative of the index of a child node  OR
//   * index of a child cell;
//   * 0 or more slots for user storage.
// A node's children are all other nodes or all cells. Node indices are stored
// as negative numbers. As an example, let ni be the index of a node. Then
// tree[ni] is the parent node index, tree[ni+1] is the number of children, and
// these are tree[ni+2 : ni+1+tree[ni+1]]. If idxs has negative numbers, then
// these are indices of more nodes; if positive, then they are indices into the
// geometric cell array geo_c2n.
//   If perform_checks, print debug messages if there are any bugs in the tree,
// and return false.

// Routines that encode the conventions above:
// Node's parent node.
template <typename IdxArrayT> KOKKOS_INLINE_FUNCTION
Int node_parent (const IdxArrayT& tree, const Int& ni)
{ return tree[ni+1]; }
// Number of children.
template <typename IdxArrayT> KOKKOS_INLINE_FUNCTION
Int node_nkids (const IdxArrayT& tree, const Int& ni)
{ return tree[ni]; }
// Whether the node's children are cells or nodes.
template <typename IdxArrayT> KOKKOS_INLINE_FUNCTION
Int node_has_cells (const IdxArrayT& tree, const Int& ni)
{ return tree[ni] > 0 && tree[ni+2] >= 0; }
// The ki'th child of node ni.
template <typename IdxArrayT> KOKKOS_INLINE_FUNCTION
Int node_kid (const IdxArrayT& tree, const Int& ni, const Int& ki)
{ const auto& kid = tree[ni+2+ki]; return kid >= 0 ? kid : -kid; }
// The user data of node ni.
template <typename IdxArrayT> KOKKOS_INLINE_FUNCTION
Int* node_slots (IdxArrayT& tree, const Int& ni)
{ return &tree[ni + 2 + tree[ni]]; }
template <typename IdxArrayT> KOKKOS_INLINE_FUNCTION
const Int* node_slots_const (IdxArrayT& tree, const Int& ni)
{ return &tree[ni + 2 + tree[ni]]; }

// Example: "(((0 1 2) (3 4)) ((5 6 7) ((8 9) ((10) (11 12)))))"
void read_tree(const char* str, std::vector<Int>& tree, const Int nslots=0);
// Read a 1D tree and output a 2D tree that is its tensor product, with 1D DOF IDs.
void read_tensor2d_tree(const char* str, std::vector<Int>& tree, const Int nslots=0);

Int get_max_tree_cell_idx(const Int* tree);

Int test_array_tree();

} // namespace tree
} // namespace slmm

#endif
