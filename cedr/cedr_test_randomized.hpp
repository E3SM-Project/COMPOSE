#ifndef INCLUDE_CEDR_TEST_RANDOMIZED_HPP
#define INCLUDE_CEDR_TEST_RANDOMIZED_HPP

#include "cedr_cdr.hpp"
#include "cedr_mpi.hpp"
#include "cedr_util.hpp"

namespace cedr {
namespace test {

class TestRandomized {
public:
  TestRandomized(const std::string& cdr_name, const mpi::Parallel::Ptr& p,
                 const Int& ncells, const bool verbose = false);

  // The subclass should call this, probably in its constructor.
  void init();

  Int run(const Int nrepeat = 1, const bool write=false);

private:
  const std::string cdr_name_;

protected:
  struct Tracer {
    typedef ProblemType PT;
    
    Int idx;
    Int problem_type;
    Int perturbation_type;
    bool no_change_should_hold, safe_should_hold, local_should_hold;
    bool write;

    std::string str() const;

    Tracer ()
      : idx(-1), problem_type(-1), perturbation_type(-1), no_change_should_hold(false),
        safe_should_hold(true), local_should_hold(true), write(false)
    {}
  };

  struct Values {
    Values (const Int ntracers, const Int ncells)
      : ncells_(ncells), v_((4*ntracers + 1)*ncells)
    {}
    Int ncells () const { return ncells_; }
    Real* rhom () { return v_.data(); }
    Real* Qm_min  (const Int& ti) { return v_.data() + ncells_*(1 + 4*ti    ); }
    Real* Qm      (const Int& ti) { return v_.data() + ncells_*(1 + 4*ti + 1); }
    Real* Qm_max  (const Int& ti) { return v_.data() + ncells_*(1 + 4*ti + 2); }
    Real* Qm_prev (const Int& ti) { return v_.data() + ncells_*(1 + 4*ti + 3); }
    const Real* rhom () const { return const_cast<Values*>(this)->rhom(); }
    const Real* Qm_min  (const Int& ti) const
    { return const_cast<Values*>(this)->Qm_min (ti); }
    const Real* Qm      (const Int& ti) const
    { return const_cast<Values*>(this)->Qm     (ti); }
    const Real* Qm_max  (const Int& ti) const
    { return const_cast<Values*>(this)->Qm_max (ti); }
    const Real* Qm_prev (const Int& ti) const
    { return const_cast<Values*>(this)->Qm_prev(ti); }
  private:
    Int ncells_;
    std::vector<Real> v_;
  };

  // For solution output, if requested.
  struct Writer {
    std::unique_ptr<FILE, cedr::util::FILECloser> fh;
    std::vector<Int> ngcis;  // Number of i'th rank's gcis_ array.
    std::vector<Long> gcis;  // Global cell indices packed by rank's gcis_ vector.
    std::vector<int> displs; // Cumsum of above.
    ~Writer();
  };

  const mpi::Parallel::Ptr p_;
  const Int ncells_;
  // Global mesh entity IDs, 1-1 with reduction array index or QLT leaf node.
  std::vector<Long> gcis_;
  std::vector<Tracer> tracers_;

  // Tell this class the CDR.
  virtual CDR& get_cdr() = 0;

  // Fill gcis_.
  virtual void init_numbering() = 0;

  // Using tracers_, the vector of Tracers, initialize the CDR's tracers.
  virtual void init_tracers() = 0;

  virtual void run_impl(const Int trial) = 0;

private:
  // For optional output.
  bool write_inited_;
  std::shared_ptr<Writer> w_; // Only on root.

  void init_tracers_vector();

  void add_const_to_Q(
    const Tracer& t, Values& v,
    // Move 0 < alpha <= 1 of the way to the QLT or safety feasibility bound.
    const Real& alpha,
    // Whether the modification should be done in a mass-conserving way.
    const bool conserve_mass,
    // Only safety problem is feasible.
    const bool safety_problem);

  void perturb_Q(const Tracer& t, Values& v);
  void init_writer();
  void gather_field(const Real* Qm_lcl, std::vector<Real>& Qm_gbl,
                    std::vector<Real>& wrk);
  void write_field(const std::string& tracer_name, const std::string& field_name,
                   const std::vector<Real>& Qm);
  void write_pre(const Tracer& t, Values& v);
  void write_post(const Tracer& t, Values& v);

  static void generate_rho(Values& v);
  static void generate_Q(const Tracer& t, Values& v);
  static void permute_Q(const Tracer& t, Values& v);
  static std::string get_tracer_name(const Tracer& t);
  static Int check(const std::string& cdr_name, const mpi::Parallel& p,
                   const std::vector<Tracer>& ts, const Values& v);
};

} // namespace test
} // namespace cedr

#endif
