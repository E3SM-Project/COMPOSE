#include "cedr_caas.hpp"
#include "cedr_util.hpp"

namespace cedr {
namespace caas {

struct OpData { int nsum, nmin, nmax; };
static OpData g_op_data;
static void all_reduce_op (Real* in, Real* inout, int* len,
                           MPI_Datatype* /*datatype*/) {
  const int n = g_op_data.nsum + g_op_data.nmin + g_op_data.nmax;
  for (int i = 0; i < *len; ++i) {
    int k = 0;
    for ( ; k < g_op_data.nsum; ++k)
      inout[k] += in[k];
    for ( ; k < g_op_data.nmin; ++k)
      inout[k] = std::min(inout[k], in[k]);
    for ( ; k < g_op_data.nmax; ++k)
      inout[k] = std::max(inout[k], in[k]);
    in += n;
    inout += n;
  }
}

template <typename ES>
CAAS<ES>::CAAS (const mpi::Parallel::Ptr& p, const Int nlclcells)
  : p_(p), nlclcells_(nlclcells), ntracers_(0), op_(all_reduce_op, true)
{
  cedr_throw_if(true, "WIP: Can't call yet.");
}

template <typename ES>
void CAAS<ES>::declare_tracer (int problem_type) {
  cedr_throw_if( ! (problem_type & ProblemType::shapepreserve) ||
                   (problem_type & ProblemType::conserve),
                 "CAAS is a WIP; only shapepreserve (=> consistent) is "
                 "supported right now.");
  ++ntracers_;
}

template <typename ES>
void CAAS<ES>::end_tracer_declarations () {
  d_ = RealList("CAAS data", nlclcells_ * (3*ntracers_ + 1));
}

template <typename ES>
int CAAS<ES>::get_problem_type (const Int& tracer_idx) const {
  return ProblemType::shapepreserve | ProblemType::consistent;
}

template <typename ES>
Int CAAS<ES>::get_num_tracers () const {
  return ntracers_;
}

template <typename ES>
void CAAS<ES>::reduce_locally () {
}

template <typename ES>
void CAAS<ES>::reduce_globally () {
  MPI_Type_contiguous(1 + 3*ntracers_, MPI_DOUBLE, &datatype_);
  MPI_Type_commit(&datatype_);
  g_op_data.nsum = 1 + ntracers_;
  g_op_data.nmin = ntracers_;
  g_op_data.nmax = ntracers_;
  int err = MPI_Allreduce(send_.data(), recv_.data(), nlclcells_, datatype_,
                          op_.get(), p_->comm());
}

template <typename ES>
void CAAS<ES>::caas () {
}

template <typename ES>
void CAAS<ES>::run () {
  reduce_locally();
  reduce_globally();
  caas();
}

} // namespace caas
} // namespace cedr
