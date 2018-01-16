For clarity, suppose your your C++ compiler is g++-4.8 in what follows. But it
can be something else.

1. Get and install the standalone Kokkos TPL:

$ git clone https://github.com/kokkos/kokkos.git
$ ./kokkos/generate_makefile.bash --with-openmp --ldflags=-fPIC --prefix=/path/to/desired/installation --compiler=g++-4.8

2. cp an existing make.inc.* file to one for your machine, say,
make.inc.mymachine. Edit it with machine-specific information. Then
    $ ln -s make.inc.machine make.inc
    $ make -j8
    $ ./siqk_runtests.py
