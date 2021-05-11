opt=-O3 -g
CXX=g++
CXXFLAGS=-DRELAX_TIME -Wstrict-overflow=0 -DSLMM_NP_GT_4 -DSLMM_NP_MAX=16

KOKKOS= # fill in with path to Kokkos installation
# Optional. Comment out if no TPL available.
#NETCDF=

SIQK=../../siqk
LINK_LAPACK_BLAS=-llapack -lblas
BUILD_DIAG=1
