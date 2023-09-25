#include "slmm_islet.hpp"
#include "slmm_util.hpp"

namespace slmm {
namespace islet {
namespace impl {

class Nodes {
  Int np, nh;
  static const bool include_bdy_val;
  std::vector<Int*> nodes;
  std::vector<Int> data, subnp;

  void set_ptrs();  

public:
  Nodes();
  Nodes(const Nodes& s);
  Nodes(const Int np_);

  void init(const Int np);
  bool init(const std::string& s);

  Int get_np () const { return np; }
  Int get_nh () const { return nh; }
  bool include_bdy () const { return include_bdy_val; }
  Int const* const* get_nodes () const { return nodes.data(); }
  const Int* get_subnp () const { return subnp.data(); }

  void set(const Int i, const std::initializer_list<Int>& il);
  void set(const Int i, const std::vector<Int>& il);
  void set(const Int ireg, const Int* const inodes, const Int isubnp);

  bool ok_to_eval() const;

  std::string string(const bool newline = true) const;
};

const bool Nodes::include_bdy_val = true;

bool operator==(const Nodes&, const Nodes&);
bool operator!=(const Nodes&, const Nodes&);

void eval(const Int& np, const Real* const xnodes,
          const Int* subnp, Int const* const* nodes,
          const Real& x, Real* const v);
void eval(const Nodes& nodes, const Real* const xnodes,
          const Real& x, Real* const v);
void eval(const Nodes& nodes, const Real& x, Real* const v);

void unittest_Nodes();

Nodes::Nodes () : np(-1) {}

Nodes::Nodes (const Nodes& s)
  : np(s.np), nh(s.nh), data(s.data), subnp(s.subnp)
{ set_ptrs(); }

Nodes::Nodes (const Int np_) { init(np_); }

void Nodes::set_ptrs () {
  nodes.resize(nh);
  for (size_t i = 0; i < nodes.size(); ++i)
    nodes[i] = data.data() + np*i;
}

void Nodes::init (const Int np_) {
  np = np_;
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

void Nodes::set (const Int ireg, const Int* const inodes, const Int isubnp) {
  assert(ireg <= static_cast<Int>(nodes.size()));
  assert(isubnp <= np);
  subnp[ireg] = isubnp;
  for (Int j = 0; j < isubnp; ++j) nodes[ireg][j] = inodes[j];
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
  if (np < 2 || np > Basis::np_max || (np > 13 && np != 16)) return false;
  Int include_bdy = read_int();
  assert(include_bdy);
  init(np);
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

void eval (const Nodes& nodes, const Real* const xnodes, const Real& x,
           Real* const v) {
  islet::eval(nodes.get_np(), xnodes,
              nodes.get_subnp(), nodes.get_nodes(),
              x, v);
}

void eval (const Nodes& nodes, const Real& x, Real* const v) {
  GLL gll;
  const Real* x_gll;
  gll.get_x(nodes.get_np(), x_gll);
  eval(nodes, x_gll, x, v);
}

} // namespace impl

void unittest_Nodes () {
  using namespace impl;
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
    GLL gll;
    const Real* x_gll;
    gll.get_x(12, x_gll);
    Nodes n1;
    n1.init("12 1 | 0 12: 0 1 2 3 4 5 6 7 8 9 10 11 | 1 9: 0 1 2 3 4 5 6 7 8 | "
            "2 10: 0 1 2 3 4 5 6 7 8 9 | 3 10: 0 1 2 3 4 5 6 7 8 9 | "
            "4 9: 1 2 3 4 5 6 7 8 9 | 5 10: 1 2 3 4 5 6 7 8 9 10");
    const Int nx = 71;
    Int ne = 0;
    for (Int ix = 0; ix <= nx; ++ix) {
      const Real x = -1 + Real(ix)/nx;
      Real v1[12], v2[12];
      evalon<12,6>(x_gll, {12, 9, 10, 10,  9, 10}, {0, 0, 0, 0, 1, 1}, x, v1);
      require(n1.ok_to_eval());
      eval(n1, x, v2);
      for (Int i = 0; i < 12; ++i) if (v1[i] != v2[i]) ++ne;
    }
    require(ne == 0);
  }
}

GllNodalFromString::GllNodalFromString () : nodes(nullptr) {}
GllNodalFromString::~GllNodalFromString () { if (nodes) delete nodes; }

bool GllNodalFromString::init (const std::string& basis) {
  nodes = new impl::Nodes();
  const bool ok = nodes->init(basis);
  if ( ! ok) return false;
  Basis::compute_weights(*this, nodes->get_np(), wts);
  return true;
}

bool GllNodalFromString::get_w (const Int& np, const Real*& wt) const {
  assert(nodes && np == nodes->get_np());
  wt = wts;
  return true;
}

bool GllNodalFromString::eval (const Int& np, const Real& x, Real* const v) const {
  assert(nodes && np == nodes->get_np());
  impl::eval(*nodes, x, v);
  return true;
}

bool GllNodalFromString::eval_derivative (const Int& np, const Real& x, Real* const v) const {
  assert(0);
  return false;  
}

FreeNodalFromString::FreeNodalFromString () : nodes(nullptr) {}
FreeNodalFromString::~FreeNodalFromString () { if (nodes) delete nodes; }

static bool read_xnodes (const Int np, const std::string& s, Real* const xnodes) {
  const auto p = s.find("x");
  if (p == std::string::npos) return false;
  std::stringstream sx(s.substr(p+1));
  for (Int i = 0; i < np; ++i) {
    if (sx.rdstate() & std::istream::eofbit) return false;
    sx >> xnodes[i];
    if (sx.rdstate() & std::istream::failbit) return false;
  }
  return true;
}

bool FreeNodalFromString::init (const std::string& basis) {
  nodes = new impl::Nodes();
  auto ok = nodes->init(basis);
  if ( ! ok) return false;
  ok = read_xnodes(nodes->get_np(), basis, xnodes);
  if ( ! ok) return false;
  Basis::compute_weights(*this, nodes->get_np(), wts);
  siqk::prarr("freenodal w", wts, nodes->get_np());
  return true;
}

bool FreeNodalFromString::get_x (const Int& np, const Real*& coord) const {
  assert(nodes && np == nodes->get_np());
  coord = xnodes;
  return true;
}

bool FreeNodalFromString::get_w (const Int& np, const Real*& wt) const {
  assert(nodes && np == nodes->get_np());
  wt = wts;
  return true;
}

bool FreeNodalFromString::eval (const Int& np, const Real& x, Real* const v) const {
  assert(nodes && np == nodes->get_np());
  impl::eval(*nodes, xnodes, x, v);
  return true;
}

bool FreeNodalFromString::eval_derivative (const Int& np, const Real& x, Real* const v) const {
  assert(0);
  return false;  
}

} // namespace islet
} // namespace slmm
