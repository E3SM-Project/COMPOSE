#include "cedr_mpi.hpp"

namespace cedr {
namespace mpi {

Parallel::Ptr make_parallel (MPI_Comm comm) {
  return std::make_shared<Parallel>(comm);
}

Int Parallel::size () const {
  int sz = 0;
  MPI_Comm_size(comm_, &sz);
  return sz;
}

Int Parallel::rank () const {
  int pid = 0;
  MPI_Comm_rank(comm_, &pid);
  return pid;
}

template <> MPI_Datatype get_type<int>() { return MPI_INT; }
template <> MPI_Datatype get_type<double>() { return MPI_DOUBLE; }
template <> MPI_Datatype get_type<long>() { return MPI_LONG_INT; }

int waitany (int count, MPI_Request* reqs, int* index, MPI_Status* stats) {
  return MPI_Waitany(count, reqs, index, stats ? stats : MPI_STATUS_IGNORE);
}

int waitall (int count, MPI_Request* reqs, MPI_Status* stats) {
  return MPI_Waitall(count, reqs, stats ? stats : MPI_STATUS_IGNORE);
}

bool all_ok (const Parallel& p, bool im_ok) {
  int ok = im_ok, msg;
  all_reduce<int>(p, &ok, &msg, 1, MPI_LAND);
  return static_cast<bool>(msg);
}

}
}
