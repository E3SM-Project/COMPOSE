#ifndef INCLUDE_MPAS_COMM_UTIL_HPP
#define INCLUDE_MPAS_COMM_UTIL_HPP

#include "mpas_mesh.hpp"
#include "mpas_comm.hpp"
#include "mpas_partition.hpp"
#include "mpas_accum.hpp"

Int get_num_own_cells(const Topo& t, const Comm& he);

void get_avg_radius_and_cell_spacing(const MpasMesh& m, Comm& he,
                                     Real& radius, Real& cell_spacing);

CellGraph::ConstPtr create_CellGraph(const std::string& mesh_fname);

// Reduce f(i, v), where i in 0:n-1 and v has size nscalar. f should set v, not
// accumulate into it. If order_invariant, wrk must have size at least
//     n*nscalar + order_invariant_accum_size_nreal*nscalar*(nthr+1).
template <typename Fn> void
repro_sum (Comm& he, const int n, const int nscalar, const Ar1& wrk,
           const Fn& f, const bool order_invariant, Real* const g)
{
  if (order_invariant) {
    Real* const v = wrk.data();
    ompparfor for (int i = 0; i < n; ++i)
      f(i, v + i*nscalar);
    const int os = n*nscalar;
    Ar1 w1(wrk.m() - os, wrk.data() + os);
    he.repro_all_reduce(Comm::SUM, n, nscalar, w1, v, g);
  } else {
    static constexpr int max_nscalar = 128; // for OpenMP reduction
    assert(nscalar < max_nscalar);
    Real l[max_nscalar] = {0};
#if defined _OPENMP && _OPENMP >= 201500
    // Need OpenMP 4.5 for array.
#   pragma omp parallel reduction (+: l)
#endif
    for (int i = 0; i < n; ++i) {
      Real v[max_nscalar];
      f(i, v);
      for (int j = 0; j < nscalar; ++j)
        l[j] += v[j];
    }
    he.all_reduce(Comm::SUM, nscalar, l, g);
  }
}

#endif
