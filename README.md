# COMPOSE
Compact Multi-moment Performance-Portable Semi-Lagrangian methods

COMPOSE provides libraries for semi-Lagrangian transport and, together or
separately, property preservation.

CEDR: Communication-Efficient Constrained Density Reconstructors.
SIQK: Sphereical Polygon Intersection and Quadrature.

First, install Kokkos:
    https://github.com/kokkos/kokkos
For example, in a typical environment using OpenMP, a simple build line is:
    ./kokkos/generate_makefile.bash --with-serial --with-openmp --prefix=/path/to/my/libs --compiler=g++
    make -j8 install

Second, configure, build, and test COMPOSE:
    cmake \
        -D Kokkos_DIR=/path/to/my/kokkos/install \
        -D CMAKE_INSTALL_PREFIX=/path/to/my/compose/install \
        /path/to/compose/repo
    make -j8
    ctest
