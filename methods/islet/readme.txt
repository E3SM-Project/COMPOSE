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
    make

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
    meam1 4.9e-15 mcV 9.6e+01 mdef 1.0e+00 w>0 1 wtr  9.44e+00 npm  7.27e-06  7.99e-06  1.62e-05 pum  5.41e-08 | np  8 subnp 6 6 7 6 offst 0 0 0 1
    meam1 4.4e-15 mcV 1.3e+02 mdef 1.0e+00 w>0 1 wtr  7.48e+00 npm  1.08e-05  9.61e-06  1.62e-05 pum  1.60e-08 | np  8 subnp 6 6 6 6 offst 0 0 0 1
    min_np  5
    meam1 4.9e-15 mcV 1.3e+02 mdef 1.0e+00 w>0 1 wtr  8.31e+00 npm  7.21e-06  8.39e-06  1.62e-05 pum  1.10e-06 | np  8 subnp 5 7 7 6 offst 0 0 0 1
    meam1 4.7e-15 mcV 2.5e+02 mdef 1.0e+00 w>0 1 wtr  1.06e+01 npm  1.46e-05  1.42e-05  2.16e-05 pum  5.78e-09 | np  8 subnp 5 5 7 6 offst 0 0 0 1
    meam1 4.2e-15 mcV 4.1e+02 mdef 1.0e+00 w>0 1 wtr  1.89e+01 npm  5.08e-05  5.51e-05  8.70e-05 pum  3.04e-09 | np  8 subnp 5 6 5 6 offst 0 0 1 1
    np  8
    min_np  6 max_np  8
    min_np  6 max_np  7
    meam1  1.55e-15 w>0 1 wtr 9.30e+00 npm 7.51e-06 8.07e-06 1.62e-05 pum  2.71e-08 | np  8 subnp 6 6 7 6 nodes | 0 1 2 3 5 6 | 0 1 2 3 4 5 | 0 1 2 3 4 5 6 | 1 2 3 4 5 6
    meam1  1.78e-15 w>0 1 wtr 9.60e+00 npm 7.71e-06 8.20e-06 1.62e-05 pum  4.29e-09 | np  8 subnp 6 6 7 6 nodes | 0 1 2 3 4 6 | 0 1 2 3 4 6 | 0 1 2 3 4 5 6 | 1 2 3 4 5 6
    meam1  1.78e-15 w>0 1 wtr 9.58e+00 npm 7.74e-06 8.21e-06 1.62e-05 pum  3.65e-08 | np  8 subnp 6 6 7 6 nodes | 0 1 2 3 4 7 | 0 1 2 3 4 6 | 0 1 2 3 4 5 6 | 1 2 3 4 5 6
    ...
    meam1  8.88e-16 w>0 1 wtr 8.74e+00 npm 9.32e-06 1.05e-05 2.13e-05 pum  1.17e-09 | np  8 subnp 6 6 7 6 nodes | 0 1 2 3 5 7 | 0 1 2 3 4 6 | 0 1 2 3 4 5 6 | 0 2 3 4 5 7
    meam1  1.55e-15 w>0 1 wtr 8.83e+00 npm 9.45e-06 1.05e-05 2.13e-05 pum  7.66e-10 | np  8 subnp 6 6 7 6 nodes | 0 1 2 3 4 7 | 0 1 2 3 4 6 | 0 1 2 3 4 6 7 | 0 2 3 4 5 7
    min_np  6 max_np  6
    meam1  1.78e-15 w>0 1 wtr 7.15e+00 npm 1.24e-05 1.07e-05 1.62e-05 pum  1.02e-09 | np  8 subnp 6 6 6 6 nodes | 0 1 2 3 4 6 | 0 1 2 3 4 6 | 0 1 2 3 4 6 | 1 2 3 4 5 6
    min_np  7 max_np  8
    min_np  7 max_np  7
    min_np  8 max_np  8

In this output, each line beginning with "meam1" corresponds to a t.p.s. basis.
"meam1" means "maximum eigenvalue amplitude minus 1", and the following value
is log10 of this quantity. It should be near machine precision. Then come a few
unused entries. Next is "w>0", which reports that all basis weights are > 0.
"wtr" is unused. "npm" lists the a_1,2,infty values for the basis. "pum" gives
the lambda_max^PUM value. After the "|" is the encoding of the basis, either
o.n.s. ("offst" is in the encoding) or general n.s.

For the methods/slmm/slmmir program, modify make.inc to point to your Kokkos
installation, then
    make
Optionally run regression tests:
    python2 slmm_runtests.py
Bash scripts in the methods/islet/figures directory call the slmmir program.

We use the language hy to create the figures. hy is a Lisp that compiles to
Python AST. We used hy 0.18.0 ('pip install hy' for the latest version) with
CPython 3.7.6 provided by Anaconda 3.

The code used to obtain performance data on Summit will be part of main E3SM
soon. The exact version used to generate the data is archived here:
    https://github.com/ambrad/E3SM/releases/tag/islet-2d-paper-summit-sl-gpu-timings
The data are here:
    https://github.com/E3SM-Project/perf-data/tree/main/nhxx-sl-summit-mar2021
