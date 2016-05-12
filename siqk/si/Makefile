opt=
CXX=g++-4.7
FC=gfortran-4.7
CXXFLAGS=$(opt) -Wall -pedantic -std=c++11 -DSIQK_FORTRAN
FFLAGS=$(opt) -ffixed-line-length-none
LDFLAGS=-lgfortran

CXXSOURCES=test.cpp fsi.cpp
F90SOURCES=testf.f90

CXXOBJECTS=$(CXXSOURCES:.cpp=.o)
F90OBJECTS=$(F90SOURCES:.f90=.o)

.SUFFIXES:
.SUFFIXES: .cpp .f .f90 .o

.cpp.o:
	$(CXX) $(CFLAGS) $(CXXFLAGS) -c $< -o $@

.f90.o:
	$(FC) $(FFLAGS) -c $< -o $@

all: test

test: $(CXXOBJECTS) $(F90OBJECTS)
	$(CXX) $(CXXOBJECTS) $(LDFLAGS) -o test.exe
	$(CXX) $(F90OBJECTS) fsi.o $(LDFLAGS) -o testf.exe

clean:
	rm -f *.o test.exe testf.exe
