#include "islet_nodalbasis.hpp"
#include "islet_npx.hpp"
#include "islet_util.hpp"

Nodes::Nodes () : np(-1), include_bdy_val(true) {}

Nodes::Nodes (const Nodes& s)
  : np(s.np), nh(s.nh), include_bdy_val(s.include_bdy_val),
    data(s.data), subnp(s.subnp)
{ set_ptrs(); }

Nodes::Nodes (const Int np_, const bool include_bdy) { init(np_, include_bdy); }

void Nodes::set_ptrs () {
  nodes.resize(nh);
  for (size_t i = 0; i < nodes.size(); ++i)
    nodes[i] = data.data() + np*i;
}

void Nodes::init (const Int np_, const bool include_bdy_) {
  np = np_;
  include_bdy_val = include_bdy_;
  nh = np/2 + (include_bdy_val ? 0 : 1);
  data.resize(np*nh);
  set_ptrs();
  subnp.resize(nh, -1);
}

void Nodes::set (const Int i, const std::initializer_list<Int>& il) {
  set(i, std::vector<Int>(il));
}

void Nodes::set (const Int i, const std::vector<Int>& il) {
  assert(i <= static_cast<Int>(nodes.size()));
  Int j = 0;
  for (const auto e : il) nodes[i][j++] = e;
  subnp[i] = j;
}

std::string Nodes::string (const bool newline) const {
  std::stringstream ss;
  ss << np << " " << int(include_bdy_val) << " | ";
  for (size_t i = 0; i < nodes.size(); ++i) {
    if (subnp[i] == -1) continue;
    ss << i << " " << subnp[i] << ":";
    for (Int j = 0; j < subnp[i]; ++j)
      ss << " " << nodes[i][j];
    if (i+1 != nodes.size()) ss << " | ";
  }
  if (newline) ss << "\n";
  return ss.str();
}

bool Nodes::init (const std::string& s) {
  std::stringstream ss(s);
  const auto read_int = [&] () -> int {
    int i;
    ss >> i;
    return i;
  };
  const auto eat_until_after_bar = [&] () {
    while ( ! ss.eof() && ss.peek() != '|') ss.get();
    ss.get();
  };
  np = read_int();
  include_bdy_val = read_int();
  init(np, include_bdy_val);
  for (Int ni = 0; ni < get_nh(); ++ni) {
    eat_until_after_bar();
    const auto ni_check = read_int();
    if (ni_check != ni) return false;
    subnp[ni] = read_int();
    if (subnp[ni] < 2 || subnp[ni] > np) return false;
    ss.get(); // colon
    for (Int i = 0; i < subnp[ni]; ++i)
      nodes[ni][i] = read_int();
  }
  return ok_to_eval();
}

bool Nodes::ok_to_eval () const {
  for (Int ni = 0; ni < get_nh(); ++ni) {
    if (subnp[ni] < 2) return false;
    for (Int i = 1; i < subnp[ni]; ++i)
      if (nodes[ni][i] <= nodes[ni][i-1])
        return false;
    Int fnd = 0;
    for (Int i = 0; i < subnp[ni]; ++i)
      if (nodes[ni][i] == ni || nodes[ni][i] == ni+1)
        ++fnd;
    if (fnd != 2) return false;
  }
  return true;
}

bool operator== (const Nodes& a, const Nodes& b) {
  if (a.get_np() != b.get_np()) return false;
  if (a.include_bdy() != b.include_bdy()) return false;
  const auto an = a.get_nodes();
  const auto bn = b.get_nodes();
  const auto as = a.get_subnp();
  const auto bs = b.get_subnp();
  for (Int ni = 0; ni < a.get_nh(); ++ni) {
    if (as[ni] != bs[ni]) return false;
    for (Int i = 0; i < as[ni]; ++i)
      if (an[ni][i] != bn[ni][i]) return false;
  }
  return true;
}

bool operator!= (const Nodes& a, const Nodes& b) { return ! (a == b); }

void eval (const Int& np, const bool bdy, const Real* const xnodes,
           const Int* subnp, Int const* const* nodes,
           const Real& x, Real* const v) {
  if (x > 0) {
    eval(np, bdy, xnodes, subnp, nodes, -x, v);
    for (int i = 0; i < np/2; ++i)
      std::swap(v[i], v[np-i-1]);
    return;
  }
  Real xsub[islet::np_max], vsub[islet::np_max];
  const Int nreg = bdy ? np-1 : np+1;
  const Int ios = bdy ? 1 : 0;
  for (Int i = 0; i < nreg; ++i) {
    if (i < np-2 && x > xnodes[i+ios]) continue;
    if (subnp[i] == np) {
      eval_lagrange_poly(xnodes, np, x, v);
    } else {
      for (Int j = 0; j < subnp[i]; ++j)
        xsub[j] = xnodes[nodes[i][j]];
      std::fill(v, v + np, 0);
      eval_lagrange_poly(xsub, subnp[i], x, vsub);
      for (Int j = 0; j < subnp[i]; ++j) {
        const auto node = nodes[i][j];
        assert(node >= 0);
        assert(node < np);
        v[node] = vsub[j];
      }
    }
    break;
  }
}

void eval (const Nodes& nodes, const Real* const xnodes, const Real& x,
           Real* const v) {
  eval(nodes.get_np(), nodes.include_bdy(), xnodes,
       nodes.get_subnp(), nodes.get_nodes(),
       x, v);
}

void eval (const Nodes& nodes, const Real& x, Real* const v) {
  eval(nodes, islet::get_x_gll_special(nodes.get_np()), x, v);
}

void unittest_Nodes () {
  Nodes n(10);
  require( ! n.ok_to_eval()); // Test all regions specified checking.
  n.set(0, {0,1,9});
  n.set(1, {1,2,9});
  require( ! n.ok_to_eval()); // Test all regions specified checking.
  n.set(2, {0,2,3,9});
  n.set(3, {0,2,3,4,9});
  n.set(4, {0,2,4,5,6});
  // Test order checking.
  require(n.ok_to_eval()); n.set(1, {2,1,9}); require( ! n.ok_to_eval()); n.set(1, {1,2,9});
  // Test interpolatory checking.
  require(n.ok_to_eval()); n.set(2, {4,5,9}); require( ! n.ok_to_eval()); n.set(2, {0,2,3,9});
  require(n.ok_to_eval());
  { Nodes n1(n); require(n1.ok_to_eval()); require(n1 == n); }
  {
    Nodes n1(n);
    n1.set(2, {1,2,3,9}  ); require(n1.ok_to_eval()); require(n1 != n);
    n1.set(2, {0,1,2,3,9}); require(n1.ok_to_eval()); require(n1 != n);
  }
  { Nodes n1; require(n1.init(n.string())); require(n1.ok_to_eval()); require(n1 == n); }
  {
    Nodes n1;
    n1.init("12 1 | 0 12: 0 1 2 3 4 5 6 7 8 9 10 11 | 1 9: 0 1 2 3 4 5 6 7 8 | "
            "2 10: 0 1 2 3 4 5 6 7 8 9 | 3 10: 0 1 2 3 4 5 6 7 8 9 | "
            "4 9: 1 2 3 4 5 6 7 8 9 | 5 10: 1 2 3 4 5 6 7 8 9 10");
    const Int nx = 71;
    Int ne = 0;
    for (Int ix = 0; ix <= nx; ++ix) {
      const Real x = -1 + Real(ix)/nx;
      Real v1[12], v2[12];
      npxstab<Real>::eval<12,6>({12, 9, 10, 10,  9, 10}, {0, 0, 0, 0, 1, 1}, x, v1);
      require(n1.ok_to_eval());
      eval(n1, x, v2);
      for (Int i = 0; i < 12; ++i)
        if (islet::reldif(v1[i], v2[i]) > 4*std::numeric_limits<Real>::epsilon()) ++ne;
    }
    require(ne == 0);
  }
}
