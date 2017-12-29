#ifndef INCLUDE_CEDR_MPI_HPP
#define INCLUDE_CEDR_MPI_HPP

#include <memory>

#include <mpi.h>

#include "cedr.hpp"

namespace cedr {
namespace mpi {

class Parallel {
  MPI_Comm comm_;
public:
  typedef std::shared_ptr<Parallel> Ptr;
  Parallel(MPI_Comm comm) : comm_(comm) {}
  MPI_Comm comm () const { return comm_; }
  Int size() const;
  Int rank() const;
  Int root () const { return 0; }
  bool amroot () const { return rank() == root(); }
};

Parallel::Ptr make_parallel(MPI_Comm comm);

template <typename T> MPI_Datatype get_type();

template <typename T>
int reduce (const Parallel& p, const T* sendbuf, T* rcvbuf, int count, MPI_Op op,
            int root) {
  MPI_Datatype dt = get_type<T>();
  return MPI_Reduce(const_cast<T*>(sendbuf), rcvbuf, count, dt, op, root, p.comm());
}

template <typename T>
int all_reduce (const Parallel& p, const T* sendbuf, T* rcvbuf, int count, MPI_Op op) {
  MPI_Datatype dt = get_type<T>();
  return MPI_Allreduce(const_cast<T*>(sendbuf), rcvbuf, count, dt, op, p.comm());
}

template <typename T>
int isend (const Parallel& p, const T* buf, int count, int dest, int tag,
           MPI_Request* ireq) {
  MPI_Datatype dt = get_type<T>();
  MPI_Request ureq;
  MPI_Request* req = ireq ? ireq : &ureq;
  int ret = MPI_Isend(const_cast<T*>(buf), count, dt, dest, tag, p.comm(), req);
  if ( ! ireq) MPI_Request_free(req);
  return ret;
}

template <typename T>
int irecv (const Parallel& p, T* buf, int count, int src, int tag, MPI_Request* ireq) {
  MPI_Datatype dt = get_type<T>();
  MPI_Request ureq;
  MPI_Request* req = ireq ? ireq : &ureq;
  int ret = MPI_Irecv(buf, count, dt, src, tag, p.comm(), req);
  if ( ! ireq) MPI_Request_free(req);
  return ret;
}

int waitany(int count, MPI_Request* reqs, int* index, MPI_Status* stats = nullptr);

int waitall(int count, MPI_Request* reqs, MPI_Status* stats = nullptr);

template<typename T>
int gather (const Parallel& p, const T* sendbuf, int sendcount,
            T* recvbuf, int recvcount, int root) {
  MPI_Datatype dt = get_type<T>();
  return MPI_Gather(sendbuf, sendcount, dt, recvbuf, recvcount, dt, root, p.comm());
}

template <typename T>
int gatherv (const Parallel& p, const T* sendbuf, int sendcount,
             T* recvbuf, const int* recvcounts, const int* displs, int root) {
  MPI_Datatype dt = get_type<T>();
  return MPI_Gatherv(sendbuf, sendcount, dt, recvbuf, recvcounts, displs, dt, root,
                     p.comm());
}

bool all_ok(const Parallel& p, bool im_ok);

}
}

#include "cedr_mpi_inl.hpp"

#endif
