#ifndef INCLUDE_SLMM_SPF_HPP
#define INCLUDE_SLMM_SPF_HPP

#include <memory>

#include "slmm_defs.hpp"

namespace slmm {
namespace spf {
namespace impl {
template <typename IdxArrayT> struct Tree;

inline Real get_xbd (const Real* xbd, const Int i, const bool xbds_scalar)
{ return xbds_scalar ? *xbd : xbd[i]; }

inline bool is_inside (const Real xi, const Real* xlo, const Real* xhi,
                       const Int i, const bool xbds_scalar) {
  return (xi > get_xbd(xlo, i, xbds_scalar) &&
          xi < get_xbd(xhi, i, xbds_scalar));
}

inline bool is_outside (const Real xi, const Real* xlo, const Real* xhi,
                       const Int i, const bool xbds_scalar) {
  return (xi < get_xbd(xlo, i, xbds_scalar) ||
          xi > get_xbd(xhi, i, xbds_scalar));
}
} // namespace

// Solve
//     min_x sum_i w(i) (x(i) - y(i))^2
//      st   a' x = b
//           xlo <= x <= xhi.
// Return 0 on success and x == y, 1 on success and x != y, -1 if infeasible, -2
// if max_its hit with no solution. See Section 3 of Bochev, Ridzal, Shashkov,
// Fast optimization-based conservative remap of scalar fields through aggregate
// mass transfer.
Int solve_1eq_bc_qp(const Int n, // length of vectors
                    const Real* w, const Real* a, const Real b,
                    // Bound constraints. if xbds_scalar, then xlo and xhi have
                    // length 1 instead of n.
                    const Real* xlo, const Real* xhi, const bool xbds_scalar,
                    const Real* y, Real* x, const Int max_its = 1000);
// Check the first-order optimality conditions. Return true if OK, false
// otherwise. If quiet, don't print anything.
bool check_1eq_bc_qp_foc(
  const char* label,  const Int n, const Real* w, const Real* a, const Real b,
  const Real* xlo, const Real* xhi, const bool xbds_scalar, const Real* y,
  const Real* x, const bool quiet=true);

// Algorithm ClipAndAssuredSum.
Int clip_and_sum(
  const Int n, const Real* w /* unused */, const Real* a, const Real b,
  const Real* xlo, const Real* xhi, const bool xbds_scalar,
  const Real* y, Real* x, const Int max_its = 1 /* unused */);

// Algorithm ClipAndAssuredGenericSum with weight vector a function of a.
Int clip_and_weighted_sum(
  const Int n, const Real* w /* unused */, const Real* a, const Real b,
  const Real* xlo, const Real* xhi, const bool xbds_scalar,
  const Real* y, Real* x, const Int max_its = 1 /* unused */);

// Returns the number of doubles needed in buf.
Int local_qlt_tensor2d_init(const Int np, std::vector<Int>& d);
Int local_qlt_tensor2d_run(
  const Int* d, Real* buf, // from init
  const Real* w, const Real* a, const Real b,
  const Real* xlo, const Real* xhi, const bool xbds_scalar,
  const Real* y, Real* x, const Int max_its = 1 /* unused */);
Int test_local_qlt_tensor2d();

// Global CAAS. Q_data are ordered (rho_mass, Q_min, Q_mass, Q_max).
void glbl_caas(const Int n, const Real* Q_data, Real* redistributed_mass,
               const Real extra_mass, Real* v);

class MassRedistributor {
public:
  typedef std::shared_ptr<MassRedistributor> Ptr;

  enum Method { qlt, caas, mn2 };

  // ne is needed just to build the tree; once the tree is built, nothing about
  // the mesh is used.
  //   redist_prop increases the amount of data to transfer by 1.5x (from 2 to 3
  // scalars per tracer per cell). It tries to spread mass proportionately to
  // what is already there; thus, e.g., if a cell currently has 0 mass, then the
  // redistributor attempts not to give it any new mass.
  MassRedistributor(const Int ne, const ARealArray& F,
                    const Method& method);

  Method get_method () const { return method_; }

  // Record during first SPF pass.
  void record(const Int ci, const Real q_min, const Real q_max, const Real* rho,
              const Real* Q);
  void record(const Int ci, const Real* q_min, const Real* q_max, const Real* rho,
              const Real* Q);

  // Call after first pass is done.
  void redistribute_mass(const Real extra_mass = 0);

  // Use these in the second pass.
  Real get_q_min (const Int ci) const { return q_lim_(ci,0); }
  Real get_q_max (const Int ci) const { return q_lim_(ci,1); }
  Real get_delta_mass(const Int ci) const;

  // For analysis.
  Real get_mass (const Int ci) const { return Q_data_(ci,2); }

  // For checking.
  Int get_np2 () const { return np2_; }

private:
  Method method_;
  const Int ncell_, np2_;
  std::shared_ptr<impl::Tree<AIdxArray> > tree_;
  AVec2s q_lim_;
  ARealArray redistributed_mass_;
  ARealArray2 Q_data_;
  ARealArray F_;
  ARealArray wrk_;

  void init_arrays();
};

Int test(const Int ne);
Int test();

} // namespace spf
} // namespace slmm

#endif
