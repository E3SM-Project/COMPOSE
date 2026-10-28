#ifndef INCLUDE_MPAS_TYPES_HPP
#define INCLUDE_MPAS_TYPES_HPP

#include <cstdio>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <vector>
#include <memory>
#include <string>

#include "mpas_macros.hpp"
#include "mpas_array.hpp"
#include "mpas_debug.hpp"

template <typename Real> struct RealTypeTraits;

template <> struct RealTypeTraits<double> {
  typedef int64_t  SameSizeInt;
  typedef uint64_t SameSizeUInt;
};

template <> struct RealTypeTraits<float> {
  typedef int32_t  SameSizeInt;
  typedef uint32_t SameSizeUInt;
};

typedef double Real;
typedef int PyInt;
typedef int Size;
typedef int Int;
typedef Real Vec[3];

typedef Array2D<Real> Ar2;
typedef Array3D<Real> Ar3;
typedef Array4D<Real> Ar4;
typedef Array2D<Int>  Ai2;
typedef Array3D<Int>  Ai3;
typedef Array2D<char> Ac2;

typedef Array2D<const Real> Acr2;
typedef Array3D<const Real> Acr3;
typedef Array4D<const Real> Acr4;
typedef Array2D<const Int>  Aci2;
typedef Array3D<const Int>  Aci3;
typedef Array2D<const char> Acc2;

typedef Array1D<Real> Ar1;
typedef Array1D<Int>  Ai1;
typedef Array1D<char> Ac1;
typedef Array1D<const Real> Acr1;
typedef Array1D<const Int>  Aci1;
typedef Array1D<const char> Acc1;

#endif
