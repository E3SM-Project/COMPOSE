include make.inc

ifdef NETCDF_FLAGS
	CXXFLAGS+=-DMPAS_HAVE_NETCDF
	LDFLAGS+=$(NETCDF_FLAGS)
endif

SOURCES=mpas_check.cpp mpas_python.cpp mpas_gallery.cpp mpas_mesh.cpp mpas_advect_ir.cpp \
	mpas_time_int.cpp mpas_util.cpp mpas_triquad.cpp mpas_siqp.cpp mpas_siqp_test.cpp \
	mpas_trajectory.cpp mpas_cedr.cpp mpas_accum.cpp mpas_partition.cpp mpas_comm_util.cpp \
	mpas_ci_util.cpp mpas_lagrangian_levels.cpp mpas_vr.cpp mpas_netcdf.cpp \
	mpas_advect_ir_test.cpp mpas_zeta.cpp mpas_order_invariant_accum.cpp mpas_comm.cpp \
	mpas_trajectory_setup.cpp mpas_trajectory_search.cpp mpas_trajectory_test.cpp

F90_SOURCES=mpas_homme_vr.F90

SERIAL_SOURCES=mpas_comm_serial.cpp mpas_debug_serial.cpp

OBJECTS=$(SOURCES:.cpp=.o) $(F90_SOURCES:.F90=.o)

SERIAL_OBJECTS=$(SERIAL_SOURCES:.cpp=.o)

MPI_SOURCES=mpas_mpi.cpp mpas_distribution.cpp mpas_mesh_bulk_data.cpp mpas_comm_mpi.cpp \
	mpas_debug_mpi.cpp

MPI_OBJECTS=$(MPI_SOURCES:.cpp=.o)

.SUFFIXES: .cpp .F90 .o

.cpp.o:
	$(LCXX) $(CFLAGS) $(CXXFLAGS) -c $< -o $@

.F90.o:
	$(F90) -g -O3 -fopenmp -fPIC -c $< -o $@

all: libmpas mpas_standalone

# Does not require MPI or NetCDF.
libmpas: LCXX=$(CXX)
libmpas: $(OBJECTS) $(SERIAL_OBJECTS)
	$(LCXX) $(OBJECTS) $(SERIAL_OBJECTS) $(LDFLAGS) -shared -o libmpas.so

# Requires MPI and NetCDF.
mpas_standalone: LCXX=$(MPICXX) -DMPAS_HAVE_MPI
mpas_standalone: $(MPI_OBJECTS) $(OBJECTS) mpas_standalone.o
	$(LCXX) $(MPI_OBJECTS) mpas_standalone.o $(OBJECTS) $(LDFLAGS) -o mpas_standalone

clean:
	rm -f $(OBJECTS) $(SERIAL_OBJECTS) $(MPI_OBJECTS) \
	libmpas.so mpas_standalone mpas_standalone.o vertremap_base.mod

# generate by running `bash make-depends.sh`
include make.depends
