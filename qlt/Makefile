include make.inc

CXXFLAGS=$(opt) -Wall -pedantic -fopenmp -std=c++11 -I$(KOKKOS)/include -DQLT_TIME
LDFLAGS=-fopenmp -L$(KOKKOS)/lib -lkokkos -ldl
LINK_LAPACK_BLAS=-llapack -lblas

SOURCES=qlt.cpp qlt_test.cpp

OBJECTS=$(SOURCES:.cpp=.o)

.cpp.o:
	$(MPICXX) $(CFLAGS) $(CXXFLAGS) -c $< -o $@

all: testqlt

testqlt: qlt_test.o qlt.o
	$(MPICXX) qlt_test.o qlt.o $(LDFLAGS) -o testqlt

clean:
	rm -f *.o testqlt

qlt.o: qlt.hpp qlt_inline.hpp qlt_kokkos.hpp
qlt_test.o: qlt.hpp
