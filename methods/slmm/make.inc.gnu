KOKKOS= # fill in with path to Kokkos installation

CXX=g++
F90=gfortran
CXXFLAGS=-O3 -g
# newer versions of kokkos:
#LDFLAGS=-L$(KOKKOS)/lib64 -lkokkoscore -ldl
# older versions of kokkos:
LDFLAGS=-L$(KOKKOS)/lib -lkokkos -ldl

# Optional. Comment out if no TPL available.
#NETCDF=

LINK_LAPACK_BLAS=-llapack -lblas

