#ifndef INCLUDE_MPAS_COMM_HPP
#define INCLUDE_MPAS_COMM_HPP

#include "mpas_types.hpp"
#include "mpas_accum.hpp"

#include <memory>

class Distribution;
class KahanBabushkaNeumaierSummer;
class LongDoubleSummer;
class MpasMesh;

struct Comm {
  typedef std::shared_ptr<Comm> Ptr;

  Comm(std::shared_ptr<const Distribution> d = nullptr,
       std::shared_ptr<const MpasMesh> m = nullptr);

  int rank() const;
  int nrank() const;
  bool amroot() const;
  bool parallel () const { return nrank() > 1; }
  int szhalo() const;

  // Distribution::cdhalo if parallel().
  const Array1D<const char>& cdhalo() const;

  std::shared_ptr<const Distribution> get_distribution() const;

  // In all of these exchange routines, the owned entities' values should be
  // filled. The exchange routine then fills all or a masked subset of the
  // haloed entities' values.

  void exchange_cell_full(
    const int nf,                   // number of fields fs
    const int* ns,                  // ns[i] is number of subfields in fs[i]
    const int nlvl,
    const Array3D<Real>* const fs); // fs[i](0:ns[i]-1, 0:nlvl-1, 0:ncell-1)

  void exchange_cell_masked(
    const Array2D<const char>& mask, // mask(0:nlvl-1, 0:nhalo-1)
    const char maskval,              // mask(li,ei) & maskval => comm fs[:](:,li,nown+ei)
    const int nf, const int* ns, const int nlvl, const Array3D<Real>* const fs);

  // On input, f(lvl,cell), whether owned or haloed, contains a 0 or non-0
  // value. On output, these values for the cell are accumulated from all ranks.
  void exchange_and_accumulate_cell(const int nlvl, const Array2D<Real>& f);
  void exchange_and_accumulate_cell(
    const int nlvl, const Array2D<OrderInvariantAccumTraits::LocalAtomic>& f);

  // Same but for edges.
  void exchange_edge_full(const int nf, const int* ns, const int nlvl,
                          const Array3D<Real>* const fs);
  void exchange_edge_masked(const Array2D<const char>& mask, const char maskval,
                            const int nf, const int* ns, const int nlvl,
                            const Array3D<Real>* const fs);

  void exchange_edge_and_cell_masked(
    const Array2D<const char>& edge_mask, const Array2D<const char>& cell_mask,
    const char maskval, const int nf, const int* ns, const int nlvl,
    const Array3D<Real>* const edge_fs, const Array3D<Real>* const cell_fs);

  enum Op { SUM, MAX, MIN };
  void all_reduce(const Op op, const Int n, const Real* l, Real* g);
  void all_reduce(const Op op, const Int n, const Int* l, Int* g);
  void all_reduce(const Op op, const Int n, const int64_t* l, int64_t* g);
  void all_reduce(const Op op, const Int n, const KahanBabushkaNeumaierSummer* l,
                  KahanBabushkaNeumaierSummer* g);
  void all_reduce(const Op op, const Int n, const LongDoubleSummer* l,
                  LongDoubleSummer* g);

  void all_reduce (const Op op, const Int n,
                   const order_invariant_accum::IntRep<double,false,false>* l,
                   order_invariant_accum::IntRep<double,false,false>* g) {
    typedef order_invariant_accum::IntRep<double,false,false>::Int IntT;
    all_reduce(op, n, reinterpret_cast<const IntT*>(l), reinterpret_cast<IntT*>(g));
  }

  void all_reduce(const Op op, const Int n,
                  const order_invariant_accum::IntRep<double,true,false>* l,
                  order_invariant_accum::IntRep<double,true,false>* g);
  
  // Summation-order-invariant reduction using OrderInvariantAccumulator.
  //   The fast index is nscalar.
  // wrk must have size at least
  //     order_invariant_accum_size_nreal*nscalar*(nthread+1).
  void repro_all_reduce(const Op op, const int nterm, const Int nscalar,
                        const Ar1& wrk, const Real* l, Real* g);

  bool test();

  static Ptr create(std::shared_ptr<const Distribution> d = nullptr,
                    std::shared_ptr<const MpasMesh> m = nullptr);

private:
  class Impl;
  std::shared_ptr<Impl> p;

  void calc_maxmag(const int nterm, const int nscalar, const Real* const v,
                   Real* maxmag);
};

// Return pointers to cell and edge GIDs. If nullptr is returned, the GIDs are
// simply 0 : Topo::ncell-1 and 0 : Topo::nedge-1. If non-null, then LID i maps
// to GID(i), where GID is the return ed array.
const Int* get_cell_gids(const Comm& c);
const Int* get_edge_gids(const Comm& c);

#endif
