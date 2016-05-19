Simple sphere interesection prototype with no-Kokkos build.

Basic build and test run:

$ g++ -std=c++11 test.cpp
$ ./a.out -n 20

For performance profiling,

$ g++ -O3 -DSIQK_TIME -std=c++11 test.cpp
$ ./a.out -n 20

You should see something like

n (-n): 20
        test_area_ot 1.276e-02 s     1.4 MB 1.228e-03 s 1.153e-02 s
true area 1.0196e+00 mesh area 1.0196e+00 relerr 3.2447e-13

The first line is the input. The second line shows total test time, memory
highwater, octree construction time, and (search, clip, and area calculation)
time. The third line shows the true overlap area, the area based on the meshes,
and the relative error. As the mesh is refined, the relative error increases
because (a) the sphere polygon area calculation is naive and (b) the edge
normals have increasing cancellation error. Each is part of the test setup and
are not used in practice.
