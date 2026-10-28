#include "mpas_comm.hpp"
#include "mpas_accum.hpp"
#include "mpas_mesh.hpp"

Comm::Comm (std::shared_ptr<const Distribution>, MpasMesh::ConstPtr) {}

int Comm::rank () const { return 0; }
int Comm::nrank () const { return 1; }
bool Comm::amroot () const { return true; }
int Comm::szhalo () const { return 0; }

const Array1D<const char>& Comm::cdhalo () const {
  static const Array1D<const char> unused;
  return unused;
}

std::shared_ptr<const Distribution>
Comm::get_distribution () const { return nullptr; }

void Comm::exchange_cell_full (const int, const int*, const int,
                               const Array3D<Real>* const) {}
void Comm::exchange_cell_masked (const Array2D<const char>& mask, const char maskval,
                                 const int nf, const int* ns, const int,
                                 const Array3D<Real>* const fs) {}
void Comm::exchange_and_accumulate_cell (const int, const Array2D<Real>&) {}
void Comm::exchange_and_accumulate_cell (
  const int, const Array2D<OrderInvariantAccumTraits::LocalAtomic>&) {}
void Comm::exchange_edge_full (const int, const int*, const int,
                               const Array3D<Real>* const) {}
void Comm::exchange_edge_masked (const Array2D<const char>&, const char,
                                 const int, const int*, const int,
                                 const Array3D<Real>* const) {}

void Comm
::exchange_edge_and_cell_masked (
  const Array2D<const char>&, const Array2D<const char>&,
  const char, const int, const int*, const int,
  const Array3D<Real>* const, const Array3D<Real>* const) {}

void Comm::all_reduce (const Op, const Int n, const Real* l, Real* g)
{ for (Int i = 0; i < n; ++i) g[i] = l[i]; }
void Comm::all_reduce (const Op, const Int n, const Int* l, Int* g)
{ for (Int i = 0; i < n; ++i) g[i] = l[i]; }
void Comm::all_reduce (const Op, const Int n, const int64_t* l, int64_t* g)
{ for (Int i = 0; i < n; ++i) g[i] = l[i]; }

void Comm
::all_reduce (const Op op, const Int n, const KahanBabushkaNeumaierSummer* l,
              KahanBabushkaNeumaierSummer* g) {
  for (Int i = 0; i < n; ++i) g[i] = l[i];
}

void Comm
::all_reduce (const Op op, const Int n, const LongDoubleSummer* l,
              LongDoubleSummer* g) {
  for (Int i = 0; i < n; ++i) g[i] = l[i];
}

void Comm
::all_reduce (const Op op, const Int n,
              const order_invariant_accum::IntRep<double,true,false>* l,
              order_invariant_accum::IntRep<double,true,false>* g) {
  for (Int i = 0; i < n; ++i) g[i] = l[i];
}

void Comm
::repro_all_reduce (const Op op, const int nterm, const Int nscalar,
                    const Ar1& wrk, const Real* l, Real* g) {
  assert(op == Comm::SUM);
  typedef order_invariant_accum::AccumType<Real,true,false>::type Accum;
  assert(wrk.m() >= 2*nscalar*Accum::nreal);
  auto* const s = reinterpret_cast<Accum*>(wrk.data());
  const auto os = nscalar*Accum::nreal;
  auto* const maxmag = wrk.data() + os;
  calc_maxmag(nterm, nscalar, l, maxmag);
  for (int j = 0; j < nscalar; ++j)
    s[j].set_maxmag(nterm*maxmag[j]);
  Ar1 w1(wrk.m() - os, wrk.data() + os);
  Accum::sum(w1, l, nterm, s, nscalar, true);
}

bool Comm::test () { return true; }

Comm::Ptr Comm::create (std::shared_ptr<const Distribution>, MpasMesh::ConstPtr) {
  return std::make_shared<Comm>();
}

const Int* get_cell_gids (const Comm& c) { return nullptr; }
const Int* get_edge_gids (const Comm& c) { return nullptr; }
