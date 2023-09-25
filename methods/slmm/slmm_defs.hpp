#ifndef INCLUDE_SLMM_DEFS_HPP
#define INCLUDE_SLMM_DEFS_HPP

#include "siqk.hpp"
#include "slmm_array.hpp"
#ifdef _OPENMP
# include <omp.h>
#endif

#define ompparfor  _Pragma("omp parallel for")
#define omppar     _Pragma("omp parallel")
#define ompfor     _Pragma("omp for")
#define ompbarrier _Pragma("omp barrier")
#define ompsingle  _Pragma("omp single")

namespace slmm {

using siqk::Int;
using siqk::Real;
typedef Int Size;

namespace ko = Kokkos;

using siqk::RawVec3s;
using siqk::RawConstVec3s;
using siqk::Vec3s;
using siqk::ConstVec3s;
using siqk::Idxs;
using siqk::ConstIdxs;
typedef ko::View<Real*[2], ko::LayoutRight, ko::MemoryTraits<ko::Unmanaged> > RawVec2s;
typedef ko::View<const Real*[2], ko::LayoutRight, ko::MemoryTraits<ko::Unmanaged> > RawConstVec2s;

// A 2D array A can be thought of as having nslices(A) rows and szslice(A)
// columns. A slice can be obtained by
//     auto ak = slice(A, k);
// We use this format for arrays of vertices and adjacency arrays, for
// example. In most or all cases, the intention is to parallelize over slices,
// so a Kokkos operator() will do work on a particular slice.
using siqk::nslices;
using siqk::szslice;
using siqk::slice;

typedef Array2Dim<Real,3> AVec3s;
typedef Array2Dim<Real,2> AVec2s;
typedef Array2D<Int> AIdxs;
typedef Array1D<Int> AIdxArray;
typedef Array1D<Real> ARealArray;
typedef Array2D<Real> ARealArray2;

inline Vec3s::HostMirror a2Vec3s (Array2Dim<Real,3>& a)
{ return Vec3s::HostMirror(a.data(), nslices(a)); }
inline Idxs::HostMirror a2Idxs (Array2D<Int>& a)
{ return Idxs::HostMirror(a.data(), nslices(a), szslice(a)); }
inline ConstVec3s::HostMirror a2ConstVec3s (const Array2Dim<Real,3>& a)
{ return Vec3s::HostMirror(const_cast<Real*>(a.data()), nslices(a)); }
inline ConstIdxs::HostMirror a2ConstIdxs (const Array2D<Int>& a)
{ return Idxs::HostMirror(const_cast<Int*>(a.data()), nslices(a), szslice(a)); }

} // namespace slmm

#endif
