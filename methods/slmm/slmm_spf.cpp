#include "slmm_gll.hpp"
#include "slmm_spf.hpp"
#include "slmm_nla.hpp"
#include "slmm_mesh.hpp"
#include "slmm_util.hpp"
#include "slmm_accum.hpp"

//extern int mods, calls;

namespace slmm {
namespace spf {
using impl::get_xbd;
using impl::is_inside;
using impl::is_outside;

static Real calc_r_tol (const Real b, const Real* a, const Real* y,
                        const Int n) {
  Real ab = std::abs(b);
  for (Int i = 0; i < n; ++i) ab = std::max(ab, std::abs(a[i]*y[i]));
  return 1e1*std::numeric_limits<Real>::epsilon()*std::abs(ab);
}

void calc_r (const Int n, const Real* w, const Real* a, const Real b,
             const Real* xlo, const Real* xhi, const bool xbds_scalar,
             const Real* y, const Real& lambda, Real* x, Real& r,
             Real& r_lambda) {
  r = 0;
  r_lambda = 0;
  for (Int i = 0; i < n; ++i) {
    const Real q = a[i]/w[i];
    const Real x_trial = y[i] + lambda*q;
    Real xtmp;
    if (x_trial < (xtmp = get_xbd(xlo, i, xbds_scalar)))
      x[i] = xtmp;
    else if (x_trial > (xtmp = get_xbd(xhi, i, xbds_scalar)))
      x[i] = xtmp;
    else {
      x[i] = x_trial;
      r_lambda += a[i]*q;
    }
    r += a[i]*x[i];
  }
  r -= b;
}

Int solve_1eq_bc_qp (const Int n, const Real* w, const Real* a, const Real b,
                     const Real* xlo, const Real* xhi, const bool xbds_scalar,
                     const Real* y, Real* x, const Int max_its) {
  const Real r_tol = calc_r_tol(b, a, y, n);

  { // Check for a quick exit.
    bool all_in = true;
    Real r = 0;
    for (Int i = 0; i < n; ++i) {
      if (is_outside(x[i], xlo, xhi, i, xbds_scalar)) {
        all_in = false;
        break;
      }
      r += a[i]*x[i];
    }
    if (all_in) {
      r -= b;
      if (std::abs(r) <= r_tol)
        return 0;
    }
  }

  { // Eval r at end points to check for feasibility, and also possibly a quick
    // exit on a common case.
    Real r = -b;
    for (Int i = 0; i < n; ++i) {
      x[i] = get_xbd(xlo, i, xbds_scalar);
      r += a[i]*x[i];
    }
    if (std::abs(r) <= r_tol) return 1;
    if (r > 0) return -1;
    r = -b;
    for (Int i = 0; i < n; ++i) {
      x[i] = get_xbd(xhi, i, xbds_scalar);
      r += a[i]*x[i];
    }
    if (std::abs(r) <= r_tol) return 1;
    if (r < 0) return -1;
  }

  { // Check for a quick exit: the bounds are so tight that the midpoint of the
    // box satisfies r_tol.
    Real r = -b;
    for (Int i = 0; i < n; ++i) {
      x[i] = 0.5*(get_xbd(xlo, i, xbds_scalar) + get_xbd(xhi, i, xbds_scalar));
      r += a[i]*x[i];
    }
    if (std::abs(r) <= r_tol) return 1;
  }

  const Real wall_dist = 1e-3;

  // Get lambda endpoints.
  Real lamlo = 0, lamhi = 0;
  for (Int i = 0; i < n; ++i) {
    const Real rq = w[i]/a[i];
    const Real lamlo_i = rq*(get_xbd(xlo, i, xbds_scalar) - y[i]);
    const Real lamhi_i = rq*(get_xbd(xhi, i, xbds_scalar) - y[i]);
    if (i == 0) {
      lamlo = lamlo_i;
      lamhi = lamhi_i;
    } else {
      lamlo = std::min(lamlo, lamlo_i);
      lamhi = std::max(lamhi, lamhi_i);
    }
  }
  const Real lamlo_feas = lamlo, lamhi_feas = lamhi;
  Real lambda = lamlo <= 0 && lamhi >= 0 ? 0 : lamlo;

  Int info = -2;

  // Bisection-safeguarded Newton iteration for r(lambda) = 0.
  bool prev_step_bisect = false;
  Int nbisect = 0;
  for (Int iteration = 0; iteration < max_its; ++iteration) {
    // Compute x, r, r_lambda.
    Real r, r_lambda;
    calc_r(n, w, a, b, xlo, xhi, xbds_scalar, y, lambda, x, r, r_lambda);
    // Is r(lambda) - b sufficiently == 0?
    if (std::abs(r) <= r_tol) {
      info = 1;
      break;
    }
    // Check if the lambda bounds are too close.
    if (nbisect > 64) {
      if (lamhi == lamhi_feas || lamlo == lamlo_feas) {
        // r isn't small enough and one lambda bound is on the feasibility
        // limit. The QP must not be feasible.
        info = -1;
        break;
      }
      info = 1;
      break;
    }
    // Adjust lambda bounds.
    if (r > 0)
      lamhi = lambda;
    else
      lamlo = lambda;
    if (r_lambda != 0) {
      // Newton step.
      lambda -= r/r_lambda;
    } else {
      // Force bisection.
      lambda = lamlo;
    }
    // Safeguard. The wall distance check assures progress, but use it only
    // every other potential bisection.
    const Real D = prev_step_bisect ? 0 : wall_dist*(lamhi - lamlo);
    if (lambda - lamlo < D || lamhi - lambda < D) {
      lambda = 0.5*(lamlo + lamhi);
      ++nbisect;
      prev_step_bisect = true;
    } else {
      prev_step_bisect = false;
    }
  }

  return info;
}

bool check_1eq_bc_qp_foc (
  const char* label, const Int n, const Real* w, const Real* a, const Real b,
  const Real* xlo, const Real* xhi, const bool xbds_scalar, const Real* y,
  const Real* x, const bool quiet)
{
  auto& os = std::cout;
  bool ok = true;
  Real xtmp;
  // Check the bound constraints.
  for (Int i = 0; i < n; ++i)
    if (x[i] < (xtmp = get_xbd(xlo, i, xbds_scalar))) {
      os << "x[" << i << "] = " << x[i] << " but x[i] - xlo[i] = " << (x[i] - xtmp) << "\n";
      ok = false;
    }
  for (Int i = 0; i < n; ++i)
    if (x[i] > (xtmp = get_xbd(xhi, i, xbds_scalar))) {
      os << "x[" << i << "] = " << x[i] << " but xhi[i] - x[i] = " << (xtmp - x[i]) << "\n";
      ok = false;
    }
  // Check the equality constraint.
  Real r = 0;
  for (Int i = 0; i < n; ++i)
    r += a[i]*x[i];
  r -= b;
  if (std::abs(r) > calc_r_tol(b, a, y, n)) {
    os << "r = " << r << "\n";
    ok = false;
  }
  // Check the gradient is 0 when projected into the constraints. Compute
  //     g = W (x - y)
  //     g_reduced = g - C ((C'C) \ (C'g))
  // where
  //     IA = I(:,A)
  //     C = [IA a],
  // and A is the active set.
  Real lambda = 0, den = 0;
  for (Int i = 0; i < n; ++i)
    if (is_inside(x[i], xlo, xhi, i, xbds_scalar)) {
      const Real gi = w[i]*(x[i] - y[i]);
      lambda += a[i]*gi;
      den += a[i]*a[i];
    }
  lambda /= den;
  Real normg = 0, normy = 0;
  for (Int i = 0; i < n; ++i) {
    normy += square(y[i]);
    if (is_inside(x[i], xlo, xhi, i, xbds_scalar))
      normg += square(w[i]*(x[i] - y[i]) - a[i]*lambda);
  }
  normy = std::sqrt(normy);
  normg = std::sqrt(normg);
  const Real gtol = 1e2*std::numeric_limits<Real>::epsilon()*normy;
  if (normg > gtol) {
    os << "norm(g) = " << normg << " gtol = " << gtol << "\n";
    ok = false;
  }
  // Check the gradient at the active boundaries.
  for (Int i = 0; i < n; ++i) {
    const bool onlo = x[i] == get_xbd(xlo, i, xbds_scalar);
    const bool onhi = onlo ? false : x[i] == get_xbd(xhi, i, xbds_scalar);
    if (onlo || onhi) {
      const Real rg = w[i]*(x[i] - y[i]) - a[i]*lambda;
      if (onlo && rg < -gtol) {
        os << "onlo but rg = " << rg << "\n";
        ok = false;
      } else if (onhi && rg > gtol) {
        os << "onhi but rg = " << rg << "\n";
        ok = false;
      }
    }
  }
  if ( ! ok) os << "label: " << label << "\n";
  return ok;
}

// The weight vectors in the following code omit a factor a[i], but the
// denominators include these. That's because in many of the calculations, a[i]
// factors out of both sides.
static Int clip_and_sum_impl (
  const Int n, const Real* a, const Real b,
  const Real* xlo, const Real* xhi, const bool xbds_scalar,
  const Real* y, Real* x, Real& m, Real* v, Real& den)
{
  assert(n <= GLL::np_max*GLL::np_max);
  Real lmass = 0, umass = 0;
  for (Int i = 0; i < n; ++i) lmass += a[i]*get_xbd(xlo, i, xbds_scalar);
  for (Int i = 0; i < n; ++i) umass += a[i]*get_xbd(xhi, i, xbds_scalar);
  Int info = 1;
  if (b < lmass || b > umass) info = -1;
  Real mass = 0;
  for (Int i = 0; i < n; ++i) mass += a[i]*y[i];
  m = b - mass;
  Real xtmp;
  bool modified = false;
  for (Int i = 0; i < n; ++i) {
    if (y[i] < (xtmp = get_xbd(xlo, i, xbds_scalar)) ||
        y[i] > (xtmp = get_xbd(xhi, i, xbds_scalar))) {
      x[i] = xtmp;
      m += (y[i] - xtmp)*a[i];
      modified = true;
    } else
      x[i] = y[i];
  }
  //calls++; if (modified) mods++;
  if (m == 0) return (modified && info == 1) ? 1 : info;
  if (m > 0)
    for (Int i = 0; i < n; ++i)
      v[i] = get_xbd(xhi, i, xbds_scalar) - x[i];
  else
    for (Int i = 0; i < n; ++i)
      v[i] = x[i] - get_xbd(xlo, i, xbds_scalar);
  den = 0;
  for (Int i = 0; i < n; ++i)
    den += v[i]*a[i];
  return info;
}

Int clip_and_sum (
  const Int n, const Real* /* unused */, const Real* a, const Real b,
  const Real* xlo, const Real* xhi, const bool xbds_scalar,
  const Real* y, Real* x, const Int /* unused */)
{
  Real m, v[GLL::np_max*GLL::np_max], fac, xtmp;
  Int info = clip_and_sum_impl(n, a, b, xlo, xhi, xbds_scalar, y, x, m, v, fac);
  if (m == 0 || fac == 0) return info;

  fac = m/fac;
  for (Int i = 0; i < n; ++i) {
    x[i] += fac*v[i];
    // Clip for numerics.
    if (x[i] < (xtmp = get_xbd(xlo, i, xbds_scalar)) ||
        x[i] > (xtmp = get_xbd(xhi, i, xbds_scalar)))
      x[i] = xtmp;
  }

  return info;
}

Int clip_and_weighted_sum (
  const Int n, const Real* /* unused */, const Real* a, const Real b,
  const Real* xlo, const Real* xhi, const bool xbds_scalar,
  const Real* y, Real* x, const Int /* unused */)
{
  Real m, v[GLL::np_max*GLL::np_max], v_den, xtmp;
  Int info = clip_and_sum_impl(
    n, a, b, xlo, xhi, xbds_scalar, y, x, m, v, v_den);
  if (m == 0 || v_den == 0) return info;

  Real w[GLL::np_max*GLL::np_max] = {0}, w_den = 0;
  Real min_mass = a[0]*y[0];
  for (Int i = 1; i < n; ++i)
    min_mass = std::min(min_mass, a[i]*y[i]);
  for (Int i = 0; i < n; ++i)
    if ((m < 0 && y[i] > get_xbd(xlo, i, xbds_scalar)) ||
        (m > 0 && y[i] < get_xbd(xhi, i, xbds_scalar))) {
      // Constant mass added to each node, if possible.
      w[i] = 1/a[i];
      w_den += w[i]*a[i];
    }

  Real alpha = 0;
  if (w_den > 0) {
    alpha = 1;
    for (Int i = 0; i < n; ++i) {
      const Real wi = w[i]/w_den, vi = v[i]/v_den;
      if (wi <= vi) continue;
      const Real
        num = get_xbd(m > 0 ? xhi : xlo, i, xbds_scalar) - x[i] - m*vi,
        den = m*(wi - vi);
      if (std::abs(num) >= std::abs(den)) continue;
      alpha = std::min(alpha, num/den);
    }
  }

  for (Int i = 0; i < n; ++i) {
    if (alpha > 0)
      x[i] += m*((1 - alpha)*(v[i]/v_den) + alpha*(w[i]/w_den));
    else
      x[i] += m*(v[i]/v_den);
    if (x[i] < (xtmp = get_xbd(xlo, i, xbds_scalar)) ||
        x[i] > (xtmp = get_xbd(xhi, i, xbds_scalar)))
      x[i] = xtmp;
  }

  return info;
}

namespace impl {

static constexpr Int max_nkids = 4;

// Tree for use in SPF.
template <typename IdxArrayT>
struct Tree {
  typedef IdxArrayT IdxArray;
  // Tree from slmm::mesh, but reindexed to use numbering rather than node
  // slots.
  IdxArrayT tree;
  // Number -> node slot map.
  IdxArrayT nbr2node;
  // Pointer to each contiguous set of numbers for a level. Level 0 is the leaf
  // level.
  IdxArrayT level_ptr;
};

// Replacements for and additions to the slmm::tree functions. These use the
// reindexed tree.
namespace reindex {
template <typename IdxArrayT> KOKKOS_INLINE_FUNCTION
Int node_level (IdxArrayT& tree, const Int& ni)
{ return tree[ni + 2 + tree[ni]]; }

template <typename IdxArrayT> KOKKOS_INLINE_FUNCTION
Int node_number (IdxArrayT& tree, const Int& ni)
{ return tree[ni + 2 + tree[ni] + 1]; }

template <typename IdxArrayT> KOKKOS_INLINE_FUNCTION
Int node_has_cells (const IdxArrayT& tree, const Int& ni)
{ return node_level(tree, ni) == 1; }

template <typename IdxArrayT> KOKKOS_INLINE_FUNCTION
Int node_kid (const IdxArrayT& tree, const Int& ni, const Int& ki)
{ return tree[ni+2+ki]; }

using tree::node_parent;
using tree::node_nkids;
using tree::node_slots;
} // namespace reindex

typedef Tree<AIdxArray> TreeHost;

namespace {
Int set_level (TreeHost::IdxArray& tree, const Int k) {
  Int& level = *tree::node_slots(tree, k);
  if (tree::node_has_cells(tree, k)) {
    level = 1;
    return 2;
  }
  level = 0;
  for (Int i = 0, nkids = tree::node_nkids(tree, k); i < nkids; ++i)
    level = std::max(level, set_level(tree, tree::node_kid(tree, k, i)));
  return level + 1;
}

void count_nodes_per_level (TreeHost::IdxArray& tree,
                            TreeHost::IdxArray& level_ptr, const Int k) {
  const Int level = *tree::node_slots(tree, k);
  ++level_ptr[level+1]; // Offset up by 1 for cumsum step later.
  const Int nkids = tree::node_nkids(tree, k);
  if (tree::node_has_cells(tree, k)) {
    assert(level == 1);
    level_ptr[1] += nkids; // Again, offset by 1.
  } else {
    for (Int i = 0; i < nkids; ++i)
      count_nodes_per_level(tree, level_ptr, tree::node_kid(tree, k, i));
  }
}

void number_nodes (TreeHost::IdxArray& tree, std::vector<Int>& numbering,
                   const Int k) {
  const Int level = tree::node_slots(tree, k)[0];
  Int& number = tree::node_slots(tree, k)[1];
  number = numbering[level];
  ++numbering[level];
  if (tree::node_has_cells(tree, k))
    return;
  for (Int i = 0, nkids = tree::node_nkids(tree, k); i < nkids; ++i)
    number_nodes(tree, numbering, tree::node_kid(tree, k, i));
}

void make_nbr2node (TreeHost::IdxArray& tree, TreeHost::IdxArray& nbr2node,
                    const Int k) {
  const Int number = tree::node_slots(tree, k)[1];
  nbr2node[number] = k;
  const Int nkids = tree::node_nkids(tree, k);
  if (tree::node_has_cells(tree, k)) {
    // Cells get their parents' slot index.
    for (Int i = 0; i < nkids; ++i)
      nbr2node[tree::node_kid(tree, k, i)] = k;
  } else {
    for (Int i = 0; i < nkids; ++i)
      make_nbr2node(tree, nbr2node, tree::node_kid(tree, k, i));
  }
}

void reindex_tree (TreeHost::IdxArray& tree, TreeHost::IdxArray& nbr2node,
                   const Int k) {
  // Parent.
  tree[k+1] = reindex::node_number(tree, tree[k+1]);
  if (tree::node_has_cells(tree, k))
    return;
  // Kids.
  for (Int i = 0, nkids = tree::node_nkids(tree, k); i < nkids; ++i) {
    const Int ki = -tree[k+2+i];
    tree[k+2+i] = reindex::node_number(tree, ki);
    reindex_tree(tree, nbr2node, ki);
  }
}

void check_reindexing (TreeHost::IdxArray& tree, TreeHost::IdxArray& nbr2node,
                       const Int k, std::vector<Int>& node_hits) {
  assert(nbr2node[reindex::node_number(tree, k)] == k);
  assert(reindex::node_level(tree, k) >= 1);
  assert(reindex::node_nkids(tree, k) <= max_nkids);
  if (k != 0) assert(node_hits[reindex::node_parent(tree, k)] == 1);
  ++node_hits[reindex::node_number(tree, k)];
  const Int nkids = reindex::node_nkids(tree, k);
  if (reindex::node_has_cells(tree, k)) {
    for (Int i = 0; i < nkids; ++i)
      ++node_hits[reindex::node_kid(tree, k, i)];
  } else {
    assert(reindex::node_level(tree, k) >= 2);
    for (Int i = 0; i < nkids; ++i)
      check_reindexing(tree, nbr2node, nbr2node[reindex::node_kid(tree, k, i)],
                       node_hits);
  }
}

Int make_tree (const Int ne, TreeHost& tree, const bool perform_checks=false) {
  Int nerr = 0;

  // Slot 1 is for level; 2 is for node number.
  mesh::make_cubedsphere_tree_over_cells(ne, tree.tree, 2);

  // Determine each node's level, where a leaf is at level 0 and a node is at 1
  // level higher than its highest-level kid.
  const Int nlevels = set_level(tree.tree, 0);
  tree.level_ptr = TreeHost::IdxArray(nlevels + 1);
  tree.level_ptr.set(0);

  // Count the number of nodes in a level.
  count_nodes_per_level(tree.tree, tree.level_ptr, 0);
  if (perform_checks && tree.level_ptr[1] != 6*ne*ne)
    ++nerr;
  // Cumsum to get pointer into each level's part of the numbering.
  for (Int i = 1; i <= nlevels; ++i)
    tree.level_ptr[i] += tree.level_ptr[i-1];

  // Number the nodes. All nodes in a level are contiguous in numbering.
  std::vector<Int> numbering(nlevels, 0);
  for (Int i = 0; i < nlevels; ++i) numbering[i] = tree.level_ptr[i];
  number_nodes(tree.tree, numbering, 0);
  if (perform_checks)
    for (Int i = 1; i < nlevels; ++i)
      if (numbering[i] != tree.level_ptr[i+1])
        ++nerr;

  // Set up map of number -> node slot.
  tree.nbr2node = TreeHost::IdxArray(tree.level_ptr[nlevels]);
  make_nbr2node(tree.tree, tree.nbr2node, 0);
  if (perform_checks) {
    Int zero_cnt = 0;
    for (Int i = 0; i < 6*ne*ne; ++i) {
      const Int k = tree.nbr2node[i];
      if ( ! tree::node_has_cells(tree.tree, k)) ++nerr;
    }
    for (Int i = 6*ne*ne; i < tree.level_ptr[nlevels]; ++i) {
      const Int k = tree.nbr2node[i];
      if (k == 0) ++zero_cnt;
      if (reindex::node_number(tree.tree, k) != i) ++nerr;
    }
    if (zero_cnt != 1) ++nerr;
  }

  // Reindex the tree so that instead of storing slots, the tree stores
  // numbers. This saves on index hops during use.
  reindex_tree(tree.tree, tree.nbr2node, 0);
  if (perform_checks) {
    std::vector<Int> node_hits(tree.nbr2node.size(), 0);
    check_reindexing(tree.tree, tree.nbr2node, 0, node_hits);
    bool all1 = true;
    for (size_t i = 0; i < node_hits.size(); ++i)
      if (node_hits[i] != 1) {
        all1 = false;
        break;
      }
    if ( ! all1) ++nerr;
  }
  
  return nerr;
}
} // namespace
} // namespace impl

namespace {
void redistribute_l2r (const impl::TreeHost& tree,
                       ARealArray2& Q_data) {
  const Int nlevels = tree.level_ptr.size() - 1;
  const Int sz = szslice(Q_data);
# pragma omp parallel
  {
    // For each non-leaf node in a level, starting with the leaves' parents:
    for (Int li = 1; li < nlevels; ++li) {
#     pragma omp for
      for (Int i = tree.level_ptr[li]; i < tree.level_ptr[li+1]; ++i) {
        // Zero this node's data.
        for (Int j = 0; j < sz; ++j) Q_data(i,j) = 0;
        // Iterate over the node's kids.
        const Int k = tree.nbr2node[i];
        const Int nkids = impl::reindex::node_nkids(tree.tree, k);
        for (Int ki = 0; ki < nkids; ++ki) {
          const Int kid_i = impl::reindex::node_kid(tree.tree, k, ki);
          // Add the kid's pair to this node's pair.
          for (Int j = 0; j < sz; ++j) Q_data(i,j) += Q_data(kid_i,j);
        }
      }
    }
  }
}

void redistribute_r2l (const impl::TreeHost& tree,
                       const ARealArray2& Q_data,
                       ARealArray& redistributed_mass,
                       const Real root_mass) {
  static Real ones[impl::max_nkids] = {0};
  if (ones[0] == 0) for (Int i = 0; i < impl::max_nkids; ++i) ones[i] = 1;

  const Int nlevels = tree.level_ptr.size() - 1;
  // Set the root's mass delta to 0 or a user-provided value.
  redistributed_mass[tree.nbr2node.size() - 1] = root_mass;
# pragma omp parallel
  {
    // For each non-leaf node in a level, starting with the root:
    for (Int li = nlevels-1; li > 0; --li) {
#     pragma omp for
      for (Int i = tree.level_ptr[li]; i < tree.level_ptr[li+1]; ++i) {
        Real y[impl::max_nkids], xlo[impl::max_nkids], xhi[impl::max_nkids],
          Qm_kids[impl::max_nkids] = {0};

        const Int k = tree.nbr2node[i];
        const Int nkids = impl::reindex::node_nkids(tree.tree, k);
        assert(nkids <= impl::max_nkids);
        const Real node_mass = redistributed_mass[i];

        Real Qm_lo = 0, Qm_hi = 0;
        { // Gather the kids' Q data.
          Real Q_total = 0;
          for (Int ki = 0; ki < nkids; ++ki) {
            const Int kid_i = impl::reindex::node_kid(tree.tree, k, ki);
            xlo[ki] = Q_data(kid_i,1);
            Qm_lo += xlo[ki];
            const Real Q_mass = Q_data(kid_i,2);
            xhi[ki] = Q_data(kid_i,3);
            Qm_hi += xhi[ki];
            const Real Q_ideal_mass =
              (Q_mass >= xlo[ki] ?
               (Q_mass <= xhi[ki] ? Q_mass : xhi[ki]) :
               xlo[ki]);
            y[ki] = Q_ideal_mass;
            Q_total += Q_ideal_mass;
          }
        }

        // Check if the node QP is feasible.
        const bool lo = Qm_lo > node_mass, hi = Qm_hi < node_mass;
        if (lo || hi) {
          // It is not, so expand the bounds. Try solving a QP that adjusts the
          // q bound that is the problem.
          Real w[impl::max_nkids], q_min[impl::max_nkids],
            q_max[impl::max_nkids], q_bnd_orig[impl::max_nkids],
            q_bnd[impl::max_nkids];
          Real Qm_bnd = 0, rhom_tot = 0;
          for (Int ki = 0; ki < nkids; ++ki) {
            const Int kid_i = impl::reindex::node_kid(tree.tree, k, ki);
            const Real rhom = Q_data(kid_i,0);
            rhom_tot += rhom;
            w[ki] = rhom;
            q_min[ki] = xlo[ki] / rhom;
            q_max[ki] = xhi[ki] / rhom;
            q_bnd_orig[ki] = lo ? q_min[ki] : q_max[ki];
            Qm_bnd += (lo ? q_min[ki] : q_max[ki])*rhom;
          }
          Real q_bnd_lim = lo ? q_min[0] : q_max[0];
          if (lo) {
            for (Int ki = 1; ki < nkids; ++ki)
              q_bnd_lim = std::min(q_bnd_lim, q_min[ki]);
            for (Int ki = 0; ki < nkids; ++ki)
              q_min[ki] = q_bnd_lim;
          } else {
            for (Int ki = 1; ki < nkids; ++ki)
              q_bnd_lim = std::max(q_bnd_lim, q_max[ki]);
            for (Int ki = 0; ki < nkids; ++ki)
              q_max[ki] = q_bnd_lim;
          }
          // Check if *this* QP is feasible.
          bool feasible; {
            Real Qm_lo = 0, Qm_hi = 0;
            for (Int ki = 0; ki < nkids; ++ki) {
              Qm_lo += q_min[ki];
              Qm_hi += q_max[ki];
            }
            feasible = Qm_lo <= Qm_bnd && Qm_bnd <= Qm_hi;
          }
          if (feasible) {
            slmm::spf::solve_1eq_bc_qp(nkids, w, w, Qm_bnd, q_min, q_max, false,
                                       q_bnd_orig, q_bnd);
          } else {
            // The QP is't feasible, so set all the bound values to the tightest
            // constant permitted.
            const Real q_bnd_const = Qm_bnd / rhom_tot;
            for (Int ki = 0; ki < nkids; ++ki)
              q_bnd[ki] = q_bnd_const;
          }
          for (Int ki = 0; ki < nkids; ++ki) {
            const Real rhom = w[ki];
            xlo[ki] = (lo ? q_bnd[ki] : q_min[ki])*rhom;
            xhi[ki] = (hi ? q_bnd[ki] : q_max[ki])*rhom;
          }
        }
        
        // Solve the QP for the masses. The (lo || hi) block above assures
        // feasibility.
        Real w[impl::max_nkids] = {0};
        // Weight favors adding extra mass to the cell with more (total, so the
        // weight is linearly invariant) mass already.
        for (Int ki = 0; ki < nkids; ++ki)
          w[ki] = 1/Q_data(ki,0);
        solve_1eq_bc_qp(nkids, w, ones, node_mass, xlo, xhi, false, y,
                        Qm_kids);
#ifdef EXPENSIVE_CHECKS
        check_1eq_bc_qp_foc("r2l", nkids, ones, ones, node_mass, xlo, xhi,
                            false, y, Qm_kids, false);
#endif

        // Tell the kids their new masses.
        for (Int ki = 0; ki < nkids; ++ki) {
          const Int kid_i = impl::reindex::node_kid(tree.tree, k, ki);
          redistributed_mass[kid_i] = Qm_kids[ki];
        }
      }
    }
  }
}

void run_mn2 (const ARealArray2& Q_data,
              ARealArray& redistributed_mass,
              const Real extra_mass, ARealArray& v) {
  const Int ncell = redistributed_mass.size();
  Real Q_mass = 0;
  for (Int ci = 0; ci < ncell; ++ci)
    Q_mass += Q_data(ci,2);
  if (v.size() != ncell)
    resize(v, 4*ncell);
  for (Int ci = 0; ci < ncell; ++ci) {
    v(ci) = 1;
    v(  ncell+ci) = Q_data(ci,1);
    v(2*ncell+ci) = Q_data(ci,2);
    v(3*ncell+ci) = Q_data(ci,3);
  }
  for (Int ci = 0; ci < ncell; ++ci)
    redistributed_mass(ci) = Q_data(ci,2);
  solve_1eq_bc_qp(ncell, &v(0), &v(0), Q_mass + extra_mass,
                  &v(ncell), &v(3*ncell), false,
                  &v(2*ncell), &redistributed_mass(0), 100);
}

void run_caas (const ARealArray2& Q_data,
               ARealArray& redistributed_mass,
               const Real extra_mass, ARealArray& v) {
  const Int n = redistributed_mass.size();
  if (v.size() != n)
    resize(v, n);
  glbl_caas(n, Q_data.data(), redistributed_mass.data(), extra_mass, v.data());
}
} // namespace
    
void glbl_caas (const Int n, const Real* Q_data_r, Real* redistributed_mass_r,
                const Real extra_mass, Real* v_r) {
  ARealArray2 Q_data(const_cast<Real*>(Q_data_r), n, 4);
  ARealArray redistributed_mass(redistributed_mass_r, n);
  ARealArray v(v_r, n);
  Real slot[slmm::accumulate_threaded_bfb_nslot];
# pragma omp parallel
  {
#   pragma omp for
    for (Int ci = 0; ci < n; ++ci) {
      if (Q_data(ci,2) < Q_data(ci,1))
        redistributed_mass(ci) = Q_data(ci,1) - Q_data(ci,2);
      else if (Q_data(ci,2) > Q_data(ci,3))
        redistributed_mass(ci) = Q_data(ci,3) - Q_data(ci,2);
      else
        redistributed_mass(ci) = 0;
    }
    const auto m = extra_mass -
      accumulate_threaded_bfb(redistributed_mass.data(), n, slot);
    if (m > 0) {
#     pragma omp for
      for (Int ci = 0; ci < n; ++ci)
        v(ci) = (Q_data(ci,2) >= Q_data(ci,3) ?
                 0 :
                 Q_data(ci,3) - (Q_data(ci,2) + redistributed_mass(ci)));
    } else {
#     pragma omp for
      for (Int ci = 0; ci < n; ++ci)
        v(ci) = (Q_data(ci,2) <= Q_data(ci,1) ?
                 0 :
                 (Q_data(ci,2) + redistributed_mass(ci)) - Q_data(ci,1));
    }
    const auto vsum = accumulate_threaded_bfb(v.data(), n, slot);
    const auto fac = m/vsum;
#   pragma omp for
    for (Int ci = 0; ci < n; ++ci)
      redistributed_mass(ci) += Q_data(ci,2) + fac*v(ci);
  }
}

MassRedistributor
::MassRedistributor (const Int ne, const ARealArray& F,
                     const Method& method)
  : method_(method), ncell_(6*ne*ne), np2_(F.size() / ncell_)
{
  if (method_ == qlt) {
    tree_ = std::make_shared<impl::TreeHost>();
    impl::make_tree(ne, *tree_);
  }
  F_ = F;
  init_arrays();
}

void MassRedistributor::init_arrays () {
  q_lim_ = AVec2s(ncell_);
  Q_data_ = ARealArray2(
    tree_ ? tree_->nbr2node.size() : ncell_, 4);
  redistributed_mass_ = ARealArray(nslices(Q_data_));
}

void MassRedistributor
::record (const Int ci, const Real q_min, const Real q_max, const Real* rho,
          const Real* Q) {
  q_lim_(ci,0) = q_min;
  q_lim_(ci,1) = q_max;
  const Real* const F = &F_[ci*np2_];
  Real rho_mass = 0, Q_mass = 0;
  const Int os = ci*np2_;
  for (Int i = 0; i < np2_; ++i)
    rho_mass += F[i] * rho[os + i];
  if (Q)
    for (Int i = 0; i < np2_; ++i)
      Q_mass += F[i] * Q[os + i];
  else
    Q_mass = rho_mass;
  Q_data_(ci,0) = rho_mass;
  Q_data_(ci,1) = q_min * rho_mass;
  Q_data_(ci,2) = Q_mass;
  Q_data_(ci,3) = q_max * rho_mass;
}

void MassRedistributor
::record (const Int ci, const Real* q_min, const Real* q_max, const Real* rho,
          const Real* Q) {
  const Real* const F = &F_[ci*np2_];
  Real rho_mass = 0, Q_mass = 0, Q_min = 0, Q_max = 0;
  const Int os = ci*np2_;
  for (Int i = 0; i < np2_; ++i) {
    const Real rhom = F[i] * rho[os + i];
    rho_mass += rhom;
    Q_min += rhom * q_min[i];
    Q_max += rhom * q_max[i];
  }
  assert(Q);
  for (Int i = 0; i < np2_; ++i)
    Q_mass += F[i] * Q[os + i];
  Q_data_(ci,0) = rho_mass;
  Q_data_(ci,1) = Q_min;
  Q_data_(ci,2) = Q_mass;
  Q_data_(ci,3) = Q_max;
  q_lim_(ci,0) = Q_min/rho_mass;
  q_lim_(ci,1) = Q_max/rho_mass;
}

void MassRedistributor::redistribute_mass (const Real extra_mass) {
  switch (method_) {
  case qlt:
    assert(tree_);
    redistribute_l2r(*tree_, Q_data_);
    redistribute_r2l(*tree_, Q_data_, redistributed_mass_,
                     Q_data_(nslices(Q_data_)-1, 2) + extra_mass);
    break;
  case caas:
    run_caas(Q_data_, redistributed_mass_, extra_mass, wrk_);
    break;
  case mn2:
    run_mn2(Q_data_, redistributed_mass_, extra_mass, wrk_);
    break;
  default:
    throw std::runtime_error("redistributed_mass: Not a valid enum.");
  }
}

Real MassRedistributor::get_delta_mass (const Int ci) const {
  return redistributed_mass_(ci) - Q_data_(ci,2);
}

static Int test_1eq_bc_qp () {
  bool verbose = true;
  Int nerr = 0;

  Int n;
  static const Int N = 16;
  Real w[N], a[N], b, xlo[N], xhi[N], y[N], x[N], al, au;
  const bool xbds_scalar = false;

  auto run = [&] () {
    for (int type = 0; type < 4; ++type) {
      if (type == 3 && n != 16) continue;
      for (Int i = 0; i < n; ++i) x[i] = 0;
      switch (type) {
      case 0: solve_1eq_bc_qp(n, w, a, b, xlo, xhi, xbds_scalar, y, x); break;
      case 1: clip_and_sum(n, nullptr, a, b, xlo, xhi, xbds_scalar, y, x); break;
      case 2: clip_and_weighted_sum(n, nullptr, a, b, xlo, xhi, xbds_scalar, y, x); break;
      case 3: {
        std::vector<Int> tree;
        Int bufsz = local_qlt_tensor2d_init(4, tree);
        std::vector<Real> buf(bufsz);
        local_qlt_tensor2d_run(tree.data(), buf.data(),
                               nullptr, a, b, xlo, xhi, xbds_scalar, y, x);
      }
      }
      if (type == 0) {
        const bool ok = check_1eq_bc_qp_foc(
          "unittest", n, w, a, b, xlo, xhi, xbds_scalar, y, x, verbose);
        if ( ! ok) ++nerr;
      }
      Real m = 0, den = 0;
      for (Int i = 0; i < n; ++i) {
        m += a[i]*x[i];
        den += std::abs(a[i]*x[i]);
        if (x[i] < xlo[i]) ++nerr;
        else if (x[i] > xhi[i]) ++nerr;
      }
      const Real rd = std::abs(b - m)/den;
      if (rd > 1e3*std::numeric_limits<Real>::epsilon()) {
        if (verbose) pr(puf(rd) pu(n) pu(b) pu(m));
        ++nerr;
      }
    }
  };
  auto gena = [&] () {
    for (Int i = 0; i < n; ++i)
      a[i] = 0.1 + urand();
  };
  auto genw = [&] () {
    for (Int i = 0; i < n; ++i)
      w[i] = 0.1 + urand();
  };
  auto genbnds = [&] () {
    al = au = 0;
    for (Int i = 0; i < n; ++i) {
      xlo[i] = urand() - 0.5;
      al += a[i]*xlo[i];
      xhi[i] = xlo[i] + urand();
      au += a[i]*xhi[i];
    }
  };
  auto genb = [&] (const bool in) {
    if (in) {
      const Real alpha = urand();
      b = alpha*al + (1 - alpha)*au;
    } else {
      if (urand() > 0.5)
        b = au + 0.01 + urand();
      else
        b = al - 0.01 - urand();
    }
  };
  auto geny = [&] (const bool in) {
    if (in) {
      for (Int i = 0; i < n; ++i) {
        const Real alpha = urand();
        y[i] = alpha*xlo[i] + (1 - alpha)*xhi[i];
      }
    } else if (urand() > 0.2) {
      for (Int i = 1; i < n; i += 2) {
        const Real alpha = urand();
        y[i] = alpha*xlo[i] + (1 - alpha)*xhi[i];
        assert(y[i] >= xlo[i] && y[i] <= xhi[i]);
      }      
      for (Int i = 0; i < n; i += 4)
        y[i] = xlo[i] - urand();
      for (Int i = 2; i < n; i += 4)
        y[i] = xhi[i] + urand();
    } else {
      for (Int i = 0; i < n; i += 2)
        y[i] = xlo[i] - urand();
      for (Int i = 1; i < n; i += 2)
        y[i] = xhi[i] + urand();
    }
  };
  auto b4y = [&] () {
    b = 0;
    for (Int i = 0; i < n; ++i)
      b += a[i]*y[i];
  };

  for (n = 2; n <= N; ++n) {
    const Int count = n == 2 ? 100 : 10;
    for (Int i = 0; i < count; ++i) {
      gena();
      genw();
      genbnds();
      genb(true);
      geny(true);
      run();
      b4y();
      run();
      genb(true);
      geny(false);
      run();
    }
  }

  return  nerr;
}

Int test (const Int ne) {
  impl::TreeHost tree;
  Int nerr = impl::make_tree(ne, tree, true);
  return nerr;
}

Int test () {
  impl::TreeHost tree;
  Int nerr = 0;
  nerr += test_1eq_bc_qp();
  nerr += test_local_qlt_tensor2d();
  return nerr;
}

} // namespace spf
} // namespace slmm
