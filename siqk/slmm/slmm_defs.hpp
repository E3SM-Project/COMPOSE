#ifndef INCLUDE_SLMM_DEFS_HPP
#define INCLUDE_SLMM_DEFS_HPP

#include "siqk.hpp"

namespace slmm {
using siqk::Int;
using siqk::Real;
typedef Int Size;

namespace ko = Kokkos;
using geometry = siqk::SphereGeometry;

using siqk::Vec3s;
using siqk::ConstVec3s;
using siqk::Idxs;
using siqk::ConstIdxs;
typedef ko::View<Int*, siqk::Layout> IdxArray;
typedef ko::View<const Int*, siqk::Layout> ConstIdxArray;
typedef ko::View<Real*, siqk::Layout> RealArray;
typedef ko::View<const Real*, siqk::Layout> ConstRealArray;
typedef ko::View<Real**, siqk::Layout> RealArray2;
typedef ko::View<const Real**, siqk::Layout> ConstRealArray2;

// A 2D array A can be thought of as having nslices(A) rows and szslice(A)
// columns. A slice can be obtained by
//     auto ak = slice(A, k);
// We use this format for arrays of vertices and adjacency arrays, for
// example. In most or all cases, the intention is to parallelize over slices,
// so a Kokkos operator() will do work on a particular slice.
using siqk::nslices;
using siqk::szslice;
using siqk::slice;
} // namespace slmm

#endif
