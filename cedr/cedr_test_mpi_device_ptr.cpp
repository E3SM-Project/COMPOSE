// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#include "cedr_mpi.hpp"

// Seg fault if MPI can't handle device pointers. Thus, we'll isolate
// this environment-related problem from true bugs inside COMPOSE.
template <typename ES>
void test_mpi_device_ptr () {
  using namespace cedr;
  static const int n = 42, tag = 24;
  Kokkos::View<int*, ES> r("r", n);
  const auto p = mpi::make_parallel(MPI_COMM_WORLD);
  mpi::Request req;
  mpi::irecv(*p, r.data(), n, p->root(), tag, &req);
  if (p->amroot()) {
    Kokkos::View<int*, ES> s("s", n);
    const auto sh = Kokkos::create_mirror_view(s);
    for (int i = 0; i < n; ++i) sh(i) = i;
    Kokkos::deep_copy(s, sh);
    for (int i = 0; i < p->size(); ++i)
      mpi::isend(*p, s.data(), n, i, tag);
  }
  mpi::waitall(1, &req);
}

int main (int argc, char** argv) {
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv); {
    test_mpi_device_ptr<Kokkos::DefaultHostExecutionSpace>();
    test_mpi_device_ptr<Kokkos::DefaultExecutionSpace>();
  } Kokkos::finalize();
  MPI_Finalize();
  // If we didn't seg fault, we passed.
  return 0;
}
