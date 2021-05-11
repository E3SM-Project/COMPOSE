#include "slmm_spf.hpp"
#include "slmm_array_tree.hpp"

namespace slmm {
namespace spf {
namespace {

const char* tree_desc (const Int np) {
  switch (np) {
  case  3: return nullptr;
  case  4: return "((0 1) (2 3))";
  case  5: return nullptr;
  case  6: return nullptr;
  case  7: return "((0 1) (2 3 4) (5 6))";
  case  8: return "(((0 1) (2 3)) ((4 5) (6 7)))";
  case  9: return nullptr;
  case 10: return "(((0 1) (2 3)) (4 5) ((6 7) (8 9)))";
  case 11: return "(((0 1) (2 3)) (4 5 6) ((7 8) (9 10)))";
  case 12: return "(((0 1) (2 3 4)) (5 6) ((7 8 9) (10 11)))";
  case 13: return "(((0 1) (2 3 4)) (5 6 7) ((8 9 10) (11 12)))";
  case 16: return "((((0 1) (2 3)) ((4 5) (6 7))) (((8 9) (10 11)) ((12 13) (14 15))))";
  default: return nullptr;
  }
}

struct ProbData {
  const Int n;
  Real* a, * xlo, * xhi, * y;
  Real b;
  Real* x;

  ProbData (const Int n_, Real* a_, Real b_, Real* xlo_, Real* xhi_,
            Real* y_, Real* x_)
    : n(n_), a(a_), xlo(xlo_), xhi(xhi_), y(y_), b(b_), x(x_)
  {}
};

struct Data {
  const Real* a, * xlo, * xhi, * y;
  const Real b;
  const bool xbds_scalar;
  Real* x;
  // tree
  const Int* t;
  Real* buf;

  Data (const Int* const t_, Real* const buf_,
        const Real* a_, const Real b_,
        const Real* xlo_, const Real* xhi_, const bool xbds_scalar_,
        const Real* y_, Real* x_)
    : a(a_), xlo(xlo_), xhi(xhi_), y(y_), b(b_), xbds_scalar(xbds_scalar_), x(x_),
      t(t_), buf(buf_)
  {}
};

void setup_buf(Int* t, const Int node, Int& bufidx) {
  *tree::node_slots(t, node) = bufidx;
  bufidx += 3;
  if (tree::node_has_cells(t, node)) return;
  const Int nkids = tree::node_nkids(t, node);
  for (Int i = 0; i < nkids; ++i)
    setup_buf(t, tree::node_kid(t, node, i), bufidx);
}

void l2r (const Data& d, const Int node) {
  const Int nkids = tree::node_nkids(d.t, node);
  Real lmass = 0, hmass = 0, mass = 0;
  const Int bufidx = *tree::node_slots_const(d.t, node);
  if (tree::node_has_cells(d.t, node)) {
    for (Int i = 0; i < nkids; ++i) {
      const Int dof = tree::node_kid(d.t, node, i);
      lmass += d.a[dof]*impl::get_xbd(d.xlo, dof, d.xbds_scalar);
      hmass += d.a[dof]*impl::get_xbd(d.xhi, dof, d.xbds_scalar);
      mass += d.a[dof]*d.y[dof];
    }
  } else {
    for (Int i = 0; i < nkids; ++i) {
      const Int knode = tree::node_kid(d.t, node, i);
      l2r(d, knode);
      const Int kbufidx = *tree::node_slots_const(d.t, knode);
      lmass += d.buf[kbufidx+0];
      hmass += d.buf[kbufidx+1];
      mass += d.buf[kbufidx+2];
    }
  }
  d.buf[bufidx+0] = lmass;
  d.buf[bufidx+1] = hmass;
  d.buf[bufidx+2] = mass;
}

void solve_node_problem(const ProbData& d, const Int solver_type) {
#if 1
  solve_1eq_bc_qp(d.n, d.a, d.a, d.b, d.xlo, d.xhi, false, d.y, d.x);
#else
  if (solver_type == 0)
    clip_and_sum(d.n, nullptr, d.a, d.b, d.xlo, d.xhi, false, d.y, d.x);
  else
    clip_and_weighted_sum(d.n, nullptr, d.a, d.b, d.xlo, d.xhi, false, d.y, d.x);
#endif
}

void r2l (const Data& d, const Int node, const Real mass) {
  const Int nkids = tree::node_nkids(d.t, node);
  const bool leaf = tree::node_has_cells(d.t, node);

  static const Int max_nkids = 16;
  assert(nkids <= max_nkids);
  Real a[max_nkids], xlo[max_nkids], xhi[max_nkids], y[max_nkids], x[max_nkids];
  ProbData pd(nkids, a, mass, xlo, xhi, y, x);

  if (leaf) {
    for (Int i = 0; i < nkids; ++i) {
      const Int dof = tree::node_kid(d.t, node, i);
      pd.a[i] = d.a[dof];
      pd.xlo[i] = impl::get_xbd(d.xlo, dof, d.xbds_scalar);
      pd.xhi[i] = impl::get_xbd(d.xhi, dof, d.xbds_scalar);
      pd.y[i] = d.y[dof];
      pd.x[i] = pd.y[i];
    }
  } else {
    for (Int i = 0; i < nkids; ++i) {
      const Int knode = tree::node_kid(d.t, node, i);
      const Int kbufidx = *tree::node_slots_const(d.t, knode);
      pd.a[i] = 1;
      pd.xlo[i] = d.buf[kbufidx+0];
      pd.xhi[i] = d.buf[kbufidx+1];
      pd.y[i] = d.buf[kbufidx+2];
      pd.x[i] = pd.y[i];
    }
  }

  solve_node_problem(pd, leaf ? 1 : 0);

  if (leaf) {
    for (Int i = 0; i < nkids; ++i) {
      const Int dof = tree::node_kid(d.t, node, i);
      d.x[dof] = pd.x[i];
    }
  } else {
    for (Int i = 0; i < nkids; ++i) {
      const Int knode = tree::node_kid(d.t, node, i);
      r2l(d, knode, pd.x[i]);
    }
  }
}

} // namespace

Int local_qlt_tensor2d_init (const Int np, std::vector<Int>& t) {
  const char* desc = tree_desc(np);
  if ( ! desc) {
    t.resize(1);
    t[0] = -np*np;
    return 1;
  }
  tree::read_tensor2d_tree(desc, t, 1);
  Int bufidx = 0;
  setup_buf(t.data(), 0, bufidx);
  return bufidx;
}

Int local_qlt_tensor2d_run (const Int* const t, Real* const buf,
                            const Real* /* unused */,
                            const Real* a, const Real b,
                            const Real* xlo, const Real* xhi, const bool xbds_scalar,
                            const Real* y, Real* x, const Int max_its) {
  if (t[0] < 0)
    return solve_1eq_bc_qp(-t[0], a, a, b, xlo, xhi, xbds_scalar, y, x);
  Data d(t, buf, a, b, xlo, xhi, xbds_scalar, y, x);
  l2r(d, 0);
  r2l(d, 0, d.b);
  return 1;
}

Int test_local_qlt_tensor2d () {
  Int nerr = 0;
  nerr += tree::test_array_tree();
  return nerr;
}

} // namespace spf
} // namespace slmm
