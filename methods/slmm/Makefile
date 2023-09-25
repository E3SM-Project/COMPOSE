include make.inc

SIQK=../../siqk
BUILD_DIAG=1

CXXFLAGS += -DRELAX_TIME -DSLMM_NP_GT_4 -DSLMM_NP_MAX=16 -fopenmp -std=c++14 -I$(SIQK) -I$(PWD) -I$(KOKKOS)/include -DSIQK_TIME -Wno-unused-function -fPIC
LDFLAGS += -fopenmp

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

# generate by running `bash make.depends`
include make.depends
