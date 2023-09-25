include make.inc

SOURCES = islet_isl.cpp islet_tables.cpp islet_npx.cpp islet_maxeigcomp.cpp islet_xnodes_metrics.cpp islet_pum.cpp islet_studymetrics.cpp islet_nodalbasis.cpp islet_np4.cpp
OBJECTS = $(SOURCES:.cpp=.o)
.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

all: libislet search np4 pum_sweep run_meam1_sweep cslunstab

libislet: $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -fopenmp -shared $(LINK_LAPACK_BLAS) -o libislet.so

search: libislet search.o
	$(CXX) search.o $(LDFLAGS) $(LINK_LAPACK_BLAS) -L. -lislet -fopenmp -o search

np4: libislet run_np4.o
	$(CXX) run_np4.o $(LDFLAGS) $(LINK_LAPACK_BLAS) -L. -lislet -fopenmp -o np4

pum_perturb_plot: libislet pum_perturb_plot.o
	$(CXX) pum_perturb_plot.o $(LDFLAGS) $(LINK_LAPACK_BLAS) -L. -lislet -fopenmp -o pum_perturb_plot

run_meam1_sweep: libislet run_meam1_sweep.o
	$(CXX) run_meam1_sweep.o $(LDFLAGS) $(LINK_LAPACK_BLAS) -L. -lislet -fopenmp -o run_meam1_sweep

pum_sweep: libislet pum_sweep.o
	$(CXX) pum_sweep.o $(LDFLAGS) $(LINK_LAPACK_BLAS) -L. -lislet -fopenmp -o pum_sweep

cslunstab: cslunstab.o
	$(CXX) cslunstab.o $(LDFLAGS) $(LINK_LAPACK_BLAS) -fopenmp -o cslunstab

pum_perturb_plot: libislet pum_perturb_plot.o
	$(CXX) pum_perturb_plot.o $(LDFLAGS) $(LINK_LAPACK_BLAS) -L. -lislet -fopenmp -o pum_perturb_plot

clean:
	rm -f *.o *.so search np4 pum_sweep pum_perturb_plot run_meam1_sweep

# generate by running `bash make.depends`
include make.depends
