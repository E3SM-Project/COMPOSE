#include "cedr_caas.hpp"
#include "cedr_util.hpp"
#include "cedr_test_randomized.hpp"

namespace cedr {
namespace caas {

template <typename ES>
CAAS<ES>::CAAS (const mpi::Parallel::Ptr& p, const Int nlclcells)
  : p_(p), nlclcells_(nlclcells), nrhomidxs_(0), need_conserve_(false)
{
  cedr_throw_if(nlclcells == 0, "CAAS does not support 0 cells on a rank.");
  tracer_decls_ = std::make_shared<std::vector<Decl> >();  
}

template <typename ES>
void CAAS<ES>::declare_tracer(int problem_type, const Int& rhomidx) {
  cedr_throw_if( ! (problem_type & ProblemType::shapepreserve),
                "CAAS is a WIP; ! shapepreserve is not supported yet.");
  cedr_throw_if(rhomidx > 0, "rhomidx > 0 is not supported yet.");
  tracer_decls_->push_back(Decl(problem_type, rhomidx));
  if (problem_type & ProblemType::conserve)
    need_conserve_ = true;
  nrhomidxs_ = std::max(nrhomidxs_, rhomidx+1);
}

template <typename ES>
void CAAS<ES>::end_tracer_declarations () {
  cedr_throw_if(tracer_decls_->size() == 0, "#tracers is 0.");
  cedr_throw_if(nrhomidxs_ == 0, "#rhomidxs is 0.");
  probs_ = IntList("CAAS probs", static_cast<Int>(tracer_decls_->size()));
  t2r_ = IntList("CAAS t2r", static_cast<Int>(tracer_decls_->size()));
  for (Int i = 0; i < probs_.extent_int(0); ++i) {
    probs_(i) = (*tracer_decls_)[i].probtype;
    t2r_(i) = (*tracer_decls_)[i].rhomidx;
  }
  tracer_decls_ = nullptr;
  // (rho, Qm, Qm_min, Qm_max, [Qm_prev])
  const Int e = need_conserve_ ? 1 : 0;
  d_ = RealList("CAAS data", nlclcells_ * ((3+e)*probs_.size() + 1));
  const auto nslots = 4*probs_.size();
  // (e'Qm_clip, e'Qm, e'Qm_min, e'Qm_max, [e'Qm_prev])
  send_ = RealList("CAAS send", nslots);
  recv_ = RealList("CAAS recv", nslots);
}

template <typename ES>
int CAAS<ES>::get_problem_type (const Int& tracer_idx) const {
  cedr_assert(tracer_idx >= 0 && tracer_idx < probs_.extent_int(0));
  return probs_[tracer_idx];
}

template <typename ES>
Int CAAS<ES>::get_num_tracers () const {
  return probs_.extent_int(0);
}

template <typename ES>
void CAAS<ES>::reduce_locally () {
  const Int nt = probs_.size();
  Int k = 0;
  Int os = nlclcells_;
  // Qm_clip
  for ( ; k < nt; ++k) {
    Real Qm_sum = 0, Qm_clip_sum = 0;
    for (Int i = 0; i < nlclcells_; ++i) {
      const Real Qm = d_(os+i);
      Qm_sum += (probs_(k) & ProblemType::conserve ?
                 d_(os + nlclcells_*3*nt + i) /* Qm_prev */ :
                 Qm);
      const Real Qm_min = d_(os + nlclcells_*  nt + i);
      const Real Qm_max = d_(os + nlclcells_*2*nt + i);
      const Real Qm_clip = cedr::impl::min(Qm_max, cedr::impl::max(Qm_min, Qm));
      Qm_clip_sum += Qm_clip;
      d_(os+i) = Qm_clip;
    }
    send_(     k) = Qm_clip_sum;
    send_(nt + k) = Qm_sum;
    os += nlclcells_;
  }
  k += nt;
  // Qm_min, Qm_max
  for ( ; k < 4*nt; ++k) {
    Real accum = 0;
    for (Int i = 0; i < nlclcells_; ++i)
      accum += d_(os+i);
    send_(k) = accum;
    os += nlclcells_;
  }
}

template <typename ES>
void CAAS<ES>::reduce_globally () {
  int err = mpi::all_reduce(*p_, send_.data(), recv_.data(), send_.size(), MPI_SUM);
  cedr_throw_if(err != MPI_SUCCESS,
                "CAAS::reduce_globally MPI_Allreduce returned " << err);
}

template <typename ES>
void CAAS<ES>::finish_locally () {
  const Int nt = probs_.size();
  Int os = nlclcells_;
  for (Int k = 0; k < nt; ++k) {
    const Real Qm_clip_sum = recv_(     k);
    const Real Qm_sum      = recv_(nt + k);
    const Real m = Qm_sum - Qm_clip_sum;
    if (m < 0) {
      const Real Qm_min_sum = recv_(2*nt + k);
      Real fac = Qm_clip_sum - Qm_min_sum;
      if (fac > 0) {
        fac = m/fac;
        for (Int i = 0; i < nlclcells_; ++i) {
          const Real Qm_min = d_(os + nlclcells_*  nt + i);
          Real& Qm = d_(os+i);
          Qm += fac*(Qm - Qm_min);
        }
      }
    } else if (m > 0) {
      const Real Qm_max_sum = recv_(3*nt + k);
      Real fac = Qm_max_sum - Qm_clip_sum;
      if (fac > 0) {
        fac = m/fac;
        for (Int i = 0; i < nlclcells_; ++i) {
          const Real Qm_max = d_(os + nlclcells_*2*nt + i);
          Real& Qm = d_(os+i);
          Qm += fac*(Qm_max - Qm);
        }
      }
    }
    os += nlclcells_;
  }
}

template <typename ES>
void CAAS<ES>::run () {
  reduce_locally();
  reduce_globally();
  finish_locally();
}

namespace test {
struct TestCAAS : public cedr::test::TestRandomized {
  typedef CAAS<Kokkos::DefaultExecutionSpace> CAAST;

  TestCAAS (const mpi::Parallel::Ptr& p, const Int& ncells, const bool verbose)
    : TestRandomized("CAAS", p, ncells, verbose),
      p_(p)
  {
    const auto np = p->size(), rank = p->rank();
    nlclcells_ = ncells / np;
    const Int todo = ncells - nlclcells_ * np;
    if (rank < todo) ++nlclcells_;
    caas_ = std::make_shared<CAAST>(p, nlclcells_);
    init();
  }

  CDR& get_cdr () override { return *caas_; }

  void init_numbering () override {
    const auto np = p_->size(), rank = p_->rank();
    Int start = 0;
    for (Int lrank = 0; lrank < rank; ++lrank)
      start += get_nllclcells(ncells_, np, lrank);
    gcis_.resize(nlclcells_);
    for (Int i = 0; i < nlclcells_; ++i)
      gcis_[i] = start + i;
  }

  void init_tracers () override {
    // CAAS doesn't yet support everything, so remove a bunch of the tracers.
    std::vector<TestRandomized::Tracer> tracers;
    Int idx = 0;
    for (auto& t : tracers_) {
      if ( ! (t.problem_type & ProblemType::shapepreserve) ||
           ! t.local_should_hold)
        continue;
      t.idx = idx++;
      tracers.push_back(t);
      caas_->declare_tracer(t.problem_type, 0);
    }
    tracers_ = tracers;
    caas_->end_tracer_declarations();
  }

  void run_impl (const Int trial) override {
    caas_->run();
  }

private:
  mpi::Parallel::Ptr p_;
  Int nlclcells_;
  CAAST::Ptr caas_;

  static Int get_nllclcells (const Int& ncells, const Int& np, const Int& rank) {
    Int nlclcells = ncells / np;
    const Int todo = ncells - nlclcells * np;
    if (rank < todo) ++nlclcells;
    return nlclcells;
  }
};

Int unittest (const mpi::Parallel::Ptr& p) {
  const auto np = p->size();
  Int nerr = 0;
  for (Int nlclcells : {1, 2, 4, 11}) {
    Long ncells = np*nlclcells;
    if (ncells > np) ncells -= np/2;
    nerr += TestCAAS(p, ncells, false).run(1, false);
  }
  return nerr;
}
} // namespace test
} // namespace caas
} // namespace cedr