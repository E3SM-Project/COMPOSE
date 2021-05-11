include make.inc

F90 = gfortran

CXXFLAGS += $(opt) -Wall -pedantic -fopenmp -std=c++11 -I$(SIQK) -I$(PWD) -I$(KOKKOS)/include -DSIQK_TIME -Wno-unused-function -fPIC
LDFLAGS += -fopenmp -L$(KOKKOS)/lib -lkokkos -ldl

OS_NAME := $(shell uname -s)

ifdef NETCDF
	CXXFLAGS+=-isystem$(NETCDF)/include -DSLMM_HAVE_NETCDF
	LDFLAGS+=-L$(NETCDF)/lib -lnetcdf_c++4 -lnetcdf
	ifeq ($(OS_NAME), Linux)
		LDFLAGS+=-Wl,-rpath=$(NETCDF)/lib
	endif
endif

ifdef BUILD_DIAG
	LDFLAGS += -lgfortran
	CXXFLAGS += -DSLMMIR_LAURITZEN_DIAG
endif

SLMM_SOURCES=slmm_mesh.cpp slmm_io.cpp slmm_time_int.cpp slmm_gallery.cpp slmm_util.cpp slmm_nla.cpp slmm_spf.cpp slmm_fit_extremum.cpp slmm_basis.cpp slmm_islet_string.cpp slmm_islet.cpp slmm_basis_reduced.cpp slmm_accum.cpp slmm_array_tree.cpp slmm_spf_lqlt.cpp slmm_vis.cpp

SLMM_OBJECTS=$(SLMM_SOURCES:.cpp=.o)

SLMMIR_SRC=slmmir_time_int.cpp slmmir_time_int_exp.cpp slmmir_remap_data.cpp slmmir_mono_data.cpp slmmir_remapper.cpp slmmir_remapper_isl.cpp slmmir_util.cpp slmmir_p_refine.cpp slmmir_lauritzen_diag.cpp slmmir_d2c.cpp slmmir_physgrid.cpp slmmir_snapshot.cpp

ifdef BUILD_DIAG
	SLMMIR_F90_SRC=slmmir_lauritzen_diag_f90.f90
else
	SLMMIR_F90_SRC=
endif

SLMMIR_OBJ=$(SLMMIR_SRC:.cpp=.o) $(SLMMIR_F90_SRC:.f90=.o)

.SUFFIXES: .cpp .f90 .o

.cpp.o:
	$(CXX) $(CFLAGS) $(CXXFLAGS) -c $< -o $@

.f90.o:
	$(F90) -c $< -o $@

all: slmm_test slmmir physgrid

slmm_test: $(SLMM_OBJECTS) slmm_test.o
	$(CXX) $(SLMM_OBJECTS) slmm_test.o $(LDFLAGS) $(LINK_LAPACK_BLAS) -o slmm_test

slmmir: $(SLMMIR_OBJ) $(SLMM_OBJECTS) slmmir.o
	$(CXX) slmmir.o $(SLMMIR_OBJ) $(SLMM_OBJECTS) $(LDFLAGS) $(LINK_LAPACK_BLAS) -o slmmir

physgrid: $(SLMMIR_OBJ) $(SLMM_OBJECTS) physgrid.o
	$(CXX) physgrid.o $(SLMMIR_OBJ) $(SLMM_OBJECTS) $(LDFLAGS) $(LINK_LAPACK_BLAS) -o physgrid

libslmm: $(SLMM_OBJECTS) slmm_c_compat.o
	$(CXX) $(SLMM_OBJECTS) slmm_c_compat.o $(LDFLAGS) -shared $(LINK_LAPACK_BLAS) -o libslmm.so

clean:
	rm -f *.o *.mod slmm_test slmmir libslmm.so *.gcov *.gcda *.gcno

slmm_test.o: slmm_defs.hpp slmm_mesh.hpp slmm_gll.hpp slmm_io.hpp slmm_time_int.hpp slmm_gallery.hpp slmm_islet.hpp $(SIQK)/siqk.hpp $(SIQK)/siqk_sqr.hpp
slmmir.o: slmm_defs.hpp slmm_util.hpp slmm_mesh.hpp slmm_gll.hpp slmm_io.hpp slmm_time_int.hpp slmm_gallery.hpp slmm_nla.hpp slmm_spf.hpp slmm_fit_extremum.hpp slmmir_time_int.hpp $(SIQK)/siqk.hpp $(SIQK)/siqk_quadrature.hpp $(SIQK)/siqk_sqr.hpp slmmir_remapper.hpp slmmir_remap_data.hpp slmmir_mono_data.hpp slmm_islet.hpp slmm_accum.hpp slmm_vis.hpp slmmir_time_int.hpp
slmmir_time_int.o: slmm_gallery.hpp slmm_time_int.hpp slmmir_p_refine.hpp slmmir_time_int.hpp
slmmir_time_int_exp.o: slmm_gallery.hpp slmm_time_int.hpp
slmm_test.o: slmm_defs.hpp slmm_util.hpp slmm_mesh.hpp slmm_gll.hpp slmm_io.hpp slmm_time_int.hpp slmm_gallery.hpp slmm_nla.hpp slmm_spf.hpp $(SIQK)/siqk.hpp $(SIQK)/siqk_sqr.hpp
slmm_mesh.o: slmm_mesh.hpp $(SIQK)/siqk.hpp slmm_array.hpp
slmm_spf.o: slmm_spf.hpp $(SIQK)/siqk.hpp
slmm_io.o: slmm_io.hpp
slmm_nla.o: slmm_nla.hpp
slmm_time_int.o: slmm_time_int.hpp
slmm_gallery.o: slmm_gallery.hpp
slmm_util.o: slmm_util.hpp
slmmir_remap_data.o: slmmir_remap_data.hpp
slmmir_p_refine.o: slmmir_p_refine.hpp
slmmir_mono_data.o: slmmir_mono_data.hpp
slmmir_remapper.o: slmmir_remapper.hpp slmmir_remap_data.hpp slmmir_mono_data.hpp slmmir_util.hpp slmmir.hpp slmmir_physgrid.hpp
slmmir_remapper_isl.o: slmmir_remapper.hpp slmmir_remap_data.hpp slmmir_mono_data.hpp slmmir_p_refine.hpp slmm_islet.hpp slmmir.hpp slmmir_physgrid.hpp slmm_accum.hpp
slmmir_util.o: slmmir_util.hpp
slmmir_p_refine.o: slmmir_p_refine.hpp slmm_islet.hpp
slmm_islet.o: slmm_islet.hpp
slmm_islet_string.o: slmm_islet.hpp
slmm_accum.o: slmm_accum.hpp
slmm_vis.o: slmm_vis.hpp
slmmir_physgrid.o: slmmir_physgrid.hpp
physgrid.o: slmmir_physgrid.hpp
slmm_basis.o: slmm_basis.hpp
slmm_basis_reduced.o: slmm_islet.hpp slmm_basis_reduced.hpp
