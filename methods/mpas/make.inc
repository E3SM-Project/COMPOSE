CXX=g++
MPICXX=mpicxx
F90=gfortran
CXXCOM=-Wall -pedantic -std=c++11 -g -fopenmp -fPIC
CXXFLAGS=$(CXXCOM) -O3 -DNDEBUG
LDFLAGS=-fopenmp -llapack -lblas -lgfortran

# If building with netcdf for output:
# NLOC=/home/ambrad/lib/netcdf4
# NETCDF_FLAGS=-L$(NLOC)/lib -lnetcdf_c++4 -lnetcdf -Wl,-rpath=$(NLOC)/lib -Wl,-rpath=$(NETCDF)/lib

