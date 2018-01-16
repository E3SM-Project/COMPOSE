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
int reduce(const Parallel& p, const T* sendbuf, T* rcvbuf, int count, MPI_Op op,
           int root);

template <typename T>
int all_reduce(const Parallel& p, const T* sendbuf, T* rcvbuf, int count, MPI_Op op);

template <typename T>
int isend(const Parallel& p, const T* buf, int count, int dest, int tag,
          MPI_Request* ireq);

template <typename T>
int irecv(const Parallel& p, T* buf, int count, int src, int tag, MPI_Request* ireq);

int waitany(int count, MPI_Request* reqs, int* index, MPI_Status* stats = nullptr);

int waitall(int count, MPI_Request* reqs, MPI_Status* stats = nullptr);

template<typename T>
int gather(const Parallel& p, const T* sendbuf, int sendcount,
           T* recvbuf, int recvcount, int root);

template <typename T>
int gatherv(const Parallel& p, const T* sendbuf, int sendcount,
            T* recvbuf, const int* recvcounts, const int* displs, int root);

bool all_ok(const Parallel& p, bool im_ok);

struct Op {
  typedef std::shared_ptr<Op> Ptr;

  Op (MPI_User_function* function, bool commute) {
    MPI_Op_create(function, static_cast<int>(commute), op_);
  }

  ~Op () { MPI_Op_free(op_); }

  MPI_Op* get () const { return op_; }

private:
  MPI_Op* op_;
};

} // namespace mpi
} // namespace cedr

#include "cedr_mpi_inl.hpp"

#endif
