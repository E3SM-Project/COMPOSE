#ifndef INCLUDE_MPAS_ORDER_INVARIANT_ACCUM_HPP
#define INCLUDE_MPAS_ORDER_INVARIANT_ACCUM_HPP

#include "mpas_order_invariant_accum_intreps.hpp"

namespace order_invariant_accum {

// Order-BFB-invariant accumulator.
template <typename IRepT>
struct Accumulator {
  typedef IRepT IRep;
  typedef typename IRep::Real Real;
  static constexpr int nreal =
    (sizeof(Accumulator) + sizeof(Real) - 1)/sizeof(Real);

  // nterm is an upper bound on the number of terms to sum.
  // max_term_mag is an upper bound on the magnitude of a term.
  Accumulator (const int nterm, const Real max_term_mag)
    : maxmag(nterm*max_term_mag), irep(0)
  { assert(nterm > 0 && max_term_mag); }

  // maxmag = nterm*max_term_mag.
  Accumulator (const Real maxmag_)
    : maxmag(maxmag_), irep(0)
  {}

  // Construct an object from its components.
  Accumulator (const Real maxmag_, const IRep irep_)
    : maxmag(maxmag_), irep(irep_)
  { assert(maxmag >= 0); }

  // sum += a.
  Accumulator& operator+= (const Real& a) {
    assert(std::abs(a) <= maxmag);
    if (a == 0) return *this; // protect against all-0 case with maxmag == 0
    const auto t = IRep(a/maxmag);
    irep += t;
    assert(std::abs(operator()()) <= maxmag);
    return *this;
  }

  // Accumulate another accumulator's value. nterm*max_term_mag must be the same
  // in both accumulators.
  Accumulator& operator+= (const Accumulator& a) {
    assert(maxmag == a.maxmag);
    irep += a.irep;
    assert(std::abs(operator()()) <= maxmag);
    return *this;
  }

  // Return sum.
  Real operator() () const { return maxmag*irep.real(); }

  // Need this for, e.g., MPI.
  Real get_maxmag () const { return maxmag; }
  IRep get_irep   () const { return irep;   }

  // Sum an array in parallel. wrk must have size at least
  //     sizeof(Accumulator::IRep)*nthr.
  // If nthr is less than the available number ofthreads, only nthr are used.
  static Accumulator sum(const Ar1& wrk, const int nterm, const Real* terms,
                         Real maxmag = -1);

  // Sum-reduce terms[i*nscalar + j] for i in [0,n), j in [0,nscalar-1].  
  //   accum[k] is *not* 0-inited for k in [0,nscalar). In fn, update accum with
  // =, not +=.
  //   fn(i, Real a[nscalar]) sets a[0:nscalar-1].
  //   wrk[nwrk] must have size at least
  //      sizeof(Accumulator)*nscalar*nthr.
  // If nthr is less than the available number of threads, only nthr are used.
  static void sum(const Ar1& wrk, const Real* terms, const Int nterm, Accumulator* accum,
                  const Int nscalar, const bool assert_on_nthr_too_small = false);

  // These are for low-level impl details.
  Accumulator () : maxmag(0), irep(0) {}
  void set_maxmag (const Real maxmag_) { maxmag = maxmag_; }
  void set_irep   (const IRep& irep_ ) { irep = irep_;     }

private:
  Real maxmag;
  IRep irep;
};

template <typename IRepT> inline
bool is_zero (const Accumulator<IRepT>& a) { return a.get_irep().is_zero(); }

// Acummulator type selector.
template <typename Real, bool extended_precision, bool atomic_accum>
struct AccumType {
  typedef Accumulator<IntRep<Real,extended_precision,atomic_accum> > type;
};

// The following utilities support templating on an accumulator, where the
// accumulator can be either an order_invariant_accum::Accumulator or a bare
// floating point number.

// Get the internal type that is used in an MPI reduction.
template <typename IRepT> IRepT
get_reduction_representation (const Accumulator<IRepT>& a) { return a.get_irep(); }
inline double get_reduction_representation (const double& a) { return a; }
inline float  get_reduction_representation (const float & a) { return a; }

// Construct an accumulator from an existing one.
template <typename IRepT> Accumulator<IRepT>
make_accum_from_existing (const Accumulator<IRepT>& a, const IRepT& r)
{ return Accumulator<IRepT>(a.get_maxmag(), r); }
inline double make_accum_from_existing (const double& a, const double& r) { return r; }
inline float  make_accum_from_existing (const float & a, const float & r) { return r; }

// Unit test Accumulator and the various internal representations.
bool unittest_Accumulator();

} // namespace order_invariant_accum

#include "mpas_order_invariant_accum_impl.hpp"

#endif
