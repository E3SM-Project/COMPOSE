#ifndef INCLUDE_ISLET_NODALBASIS_HPP
#define INCLUDE_ISLET_NODALBASIS_HPP

#include "islet_types.hpp"
#include "islet_util.hpp"

class Nodes {
  Int np, nh;
  bool include_bdy_val;
  std::vector<Int*> nodes;
  std::vector<Int> data, subnp;

  void set_ptrs();  

public:
  typedef std::shared_ptr<Nodes> Ptr;

  Nodes();
  Nodes(const Nodes& s);
  Nodes(const Int np_, const bool include_bdy = true);

  void init(const Int np, const bool include_bdy);
  bool init(const std::string& s);

  Int get_np () const { return np; }
  Int get_nh () const { return nh; }
  bool include_bdy () const { return include_bdy_val; }
  Int const* const* get_nodes () const { return nodes.data(); }
  const Int* get_subnp () const { return subnp.data(); }

  void set(const Int i, const std::initializer_list<Int>& il);
  void set(const Int i, const std::vector<Int>& il);

  template <typename IntT>
  void set(const Int ireg, const IntT* const inodes, const Int isubnp) {
    assert(ireg <= static_cast<Int>(nodes.size()));
    assert(isubnp <= np);
    subnp[ireg] = isubnp;
    for (Int j = 0; j < isubnp; ++j) nodes[ireg][j] = inodes[j];
  }

  bool ok_to_eval() const;

  std::string string(const bool newline = true) const;
};

bool operator==(const Nodes&, const Nodes&);
bool operator!=(const Nodes&, const Nodes&);

void eval(const Int& np, const bool bdy, const Real* const xnodes,
          const Int* subnp, Int const* const* nodes,
          const Real& x, Real* const v);
void eval(const Nodes& nodes, const Real* const xnodes,
          const Real& x, Real* const v);
void eval(const Nodes& nodes, const Real& x, Real* const v);

void unittest_Nodes();

#endif
