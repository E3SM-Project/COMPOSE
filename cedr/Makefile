include make.inc

CXXFLAGS=$(opt) -Wall -pedantic -fopenmp -std=c++11 -I$(KOKKOS)/include -DQLT_TIME
LDFLAGS=-fopenmp -L$(KOKKOS)/lib -lkokkos -ldl
LINK_LAPACK_BLAS=-llapack -lblas

SOURCES=cedr_mpi.cpp cedr_util.cpp cedr_qlt.cpp cedr_caas.cpp cedr_local.cpp \
	cedr_test.cpp cedr_test_randomized.cpp cedr_test_1d_transport.cpp

OBJECTS=$(SOURCES:.cpp=.o)

.cpp.o:
	$(MPICXX) $(CFLAGS) $(CXXFLAGS) -c $< -o $@

all: testcedr

testcedr: cedr_test.o $(OBJECTS)
	$(MPICXX) $(OBJECTS) $(LDFLAGS) -o testcedr

clean:
	rm -f *.o testcedr

cedr_qlt.o: cedr_mpi.hpp cedr_mpi_inl.hpp cedr_local.hpp cedr_local_inl.hpp \
	cedr_qlt.hpp cedr_qlt_inl.hpp cedr_kokkos.hpp cedr_util.hpp \
	cedr_test_randomized.hpp
cedr_caas.o: cedr_mpi.hpp cedr_mpi_inl.hpp cedr_local.hpp cedr_local_inl.hpp \
	cedr_caas.hpp cedr_caas_inl.hpp cedr_kokkos.hpp cedr_util.hpp \
	cedr_test_randomized.hpp
cedr_test.o: cedr_qlt.hpp cedr_util.hpp
cedr_test_1d_transport.o: cedr_qlt.hpp cedr_util.hpp
cedr_local.o: cedr_local_inl.hpp
cedr_test_randomized.o: cedr_qlt.hpp cedr_util.hpp cedr_test_randomized.hpp
