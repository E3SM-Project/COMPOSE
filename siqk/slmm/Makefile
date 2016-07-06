opt=
CXX=g++-4.7

KOKKOS=/home/ambradl/lib/kokkos/cpu
SIQK=..
LINK_LAPACK_BLAS=-llapack -lblas
# Optional. Comment out if no TPL available.
NETCDF=/home/ambradl/lib/netcdf

# Should not have to change the rest.

CXXFLAGS=$(opt) -Wall -pedantic -fopenmp -std=c++11 -I$(SIQK) -I$(KOKKOS)/include -DSIQK_TIME
LDFLAGS=-fopenmp -L$(KOKKOS)/lib -lkokkos -ldl

ifdef NETCDF
	CXXFLAGS+=-I$(NETCDF)/include -DSLMM_HAVE_NETCDF
	LDFLAGS+=-L$(NETCDF)/lib -lnetcdf_c++ -lnetcdf -Wl,-rpath=$(NETCDF)/lib
endif

SOURCES=slmm_mesh.cpp slmm_io.cpp slmm_time_int.cpp slmm_gallery.cpp slmm_util.cpp

OBJECTS=$(SOURCES:.cpp=.o)

.cpp.o:
	$(CXX) $(CFLAGS) $(CXXFLAGS) -c $< -o $@

all: slmm_test slmmir

slmm_test: $(OBJECTS) slmm_test.o
	$(CXX) $(OBJECTS) slmm_test.o $(LDFLAGS) -o slmm_test

slmmir: $(OBJECTS) slmmir.o
	$(CXX) $(OBJECTS) slmmir.o $(LDFLAGS) $(LINK_LAPACK_BLAS) -o slmmir

clean:
	rm -f *.o slmm_test slmmir

$(SIQK)/siqk.hpp: $(SIQK)/siqk_intersect.hpp $(SIQK)/siqk_geometry.hpp $(SIQK)/siqk_sqr.hpp $(SIQK)/siqk_search.hpp $(SIQK)/siqk_quadrature.hpp
slmm_test.o: slmm_defs.hpp slmm_mesh.hpp slmm_gll.hpp slmm_io.hpp slmm_time_int.hpp slmm_gallery.hpp $(SIQK)/siqk.hpp
slmmir.o: slmm_defs.hpp slmm_util.hpp slmm_mesh.hpp slmm_gll.hpp slmm_io.hpp slmm_time_int.hpp slmm_gallery.hpp $(SIQK)/siqk.hpp
slmm_mesh.o: slmm_mesh.hpp $(SIQK)/siqk.hpp
slmm_io.o: slmm_io.hpp
slmm_time_int.o: slmm_time_int.hpp
slmm_gallery.o: slmm_gallery.hpp
slmm_util.o: slmm_util.hpp
