This directory and the directory methods/slmm contain the code used to generate
the results in the Islet 2D paper, except the E3SM code for the GPU performance
results.

The directory methods/islet/figures contains the scripts used to generate data
and figures. The file "figures/figs.tex" contains the latex for the figures.
Comments before each figure explain how to generate the data and then the figure
from these data. Bash and hy scripts are those in the "figures" directory.

Programs need BLAS, LAPACK, and for slmmir, Kokkos
(https://github.com/kokkos/kokkos). NetCDF is optional and was not used for the
results in this paper. We used Kokkos version 3.1
(https://github.com/kokkos/kokkos/tree/3.3.01) in our build.

For the methods/islet programs, on a standard Linux system with GNU compiler suite,
    ln -s make.inc.gnu make.inc
    make -j8

The program "cslunstab" demonstrates the unstable classical cubic interpolation
semi-Lagrangian instances. Running it should produce no output, as in this case
all assertions pass. The program is self-contained and is meant to be read. See
the top of cslunstab.cpp for instructions.

The program "search" is used to find the Islet bases. Run as follows, in this
example for np = 8:
    OMP_NUM_THREADS=48 KMP_AFFINITY=balanced ./search findnodal_given_bestosn 8
This produces output of the following form:

    np  8
    min_np  8
    min_np  7
    min_np  6
    meam1 4.7e-15 mcV 1.3e+02 mdef 1.0e+00 w>0 1 wtr  7.48e+00 npm  1.08e-05  9.61e-06  1.62e-05 pum  4.38e-05 | np  8 subnp 6 6 6 6 offst 0 0 0 1
    meam1 4.0e-15 mcV 9.6e+01 mdef 1.0e+00 w>0 1 wtr  9.44e+00 npm  7.27e-06  7.99e-06  1.62e-05 pum  1.82e-07 | np  8 subnp 6 6 7 6 offst 0 0 0 1
    count 42592
    np  8
    NsbSearchAtom::eval 4096/42592 (  9.6%)
    meam1  4.88e-15 w>0 1 wtr 5.93e+00 npm 1.89e-05 1.60e-05 2.13e-05 pum  1.58e-07 | np  8 subnp 6 6 6 6 nodes | 0 1 2 4 5 7 | 0 1 2 3 5 6 | 0 1 2 3 5 6 | 0 2 3 4 5 7
    meam1  4.66e-15 w>0 1 wtr 6.94e+00 npm 1.34e-05 1.13e-05 1.62e-05 pum  1.39e-07 | np  8 subnp 6 6 6 6 nodes | 0 1 2 3 5 7 | 0 1 2 3 4 7 | 0 1 2 3 4 7 | 1 2 3 4 5 6
    ...
    meam1  4.44e-15 w>0 1 wtr 8.55e+00 npm 6.07e-06 7.62e-06 1.62e-05 pum  1.12e-07 | np  8 subnp 6 7 7 6 nodes | 0 1 2 3 4 5 | 0 1 2 3 4 6 7 | 0 1 2 3 4 5 6 | 1 2 3 4 5 6

In this output, each line beginning with "meam1" corresponds to a t.p.s. basis.
"meam1" means "maximum eigenvalue amplitude minus 1", and the following value
is log10 of this quantity. It should be near machine precision. Then come a few
unused entries. Next is "w>0", which reports that all basis weights are > 0.
"wtr" is unused. "npm" lists the a_1,2,infty values for the basis. "pum" gives
the lambda_max^PUM value. After the "|" is the encoding of the basis, either
o.n.s. ("offst" is in the encoding) or general n.s.

For the methods/slmm/slmmir program, modify make.inc.gnu to point to your Kokkos
installation, then
    ln -s make.inc.gnu make.inc
    make -j16
Optionally run regression tests:
    OMP_NUM_THREADS=16 KMP_AFFINITY=balanced python slmm_runtests.py
Bash scripts in the methods/islet/figures directory call the slmmir program.

We use the language hy to create the figures. hy is a Lisp that compiles to
Python AST. We used hy 0.18.0 ('pip install hy' for the latest version) with
CPython 3.7.6 provided by Anaconda 3.

The code used to obtain performance data on Summit is part of the main E3SM
repo. The exact version used to generate the data is archived here:
    https://github.com/ambrad/E3SM/releases/tag/islet-2d-paper-summit-sl-gpu-timings
The data are here:
    https://github.com/E3SM-Project/perf-data/tree/main/nhxx-sl-summit-mar2021
