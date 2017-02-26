include make.inc

CXXFLAGS=$(opt) -Wall -pedantic -fopenmp -std=c++11 -I$(KOKKOS)/include -DSIQK_TIME -Wno-unused-function
LDFLAGS=-fopenmp -L$(KOKKOS)/lib -lkokkos -ldl

OBJECTS=$(SOURCES:.cpp=.o)

.cpp.o:
	$(CXX) $(CFLAGS) $(CXXFLAGS) -c $< -o $@

all: siqk_test

siqk_test: $(OBJECTS) siqk_test.o
	$(CXX) $(OBJECTS) siqk_test.o $(LDFLAGS) -o siqk_test

clean:
	rm -f *.o siqk_test
