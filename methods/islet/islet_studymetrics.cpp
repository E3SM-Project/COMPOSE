#include "islet_pum.hpp"
#include "islet_xnodes_metrics.hpp"

namespace {

Nodes make_offset_nodal (int np, int n, const int* subnp, const int* offst) {
  Nodes nodes(np);
  const auto nh = nodes.get_nh();
  std::vector<Int> ns;
  ns.reserve(np);
  for (Int ireg = 0; ireg < nh; ++ireg) {
    if (ireg < n) {
      ns.resize(subnp[ireg]);
      for (Int i = 0; i < subnp[ireg]; ++i)
        ns[i] = offst[ireg] + i;
    } else {
      ns.resize(np);
      for (Int i = 0; i < np; ++i) ns[i] = i;
    }
    nodes.set(ireg, ns);
  }
  return nodes;
}

bool read_xnodes (const Int np, const std::string& s, Real* const xnodes) {
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

class Basis : public UserInterpMethod {
  Nodes nodes;
  bool ok, free_nodal;
  Real xnodes[islet::np_max];

public:
  Basis (const std::string& basis) {
    ok = nodes.init(basis);
    if ( ! ok) {
      printf("Invalid basis string: %s\n", basis.c_str());
      return;
    }
    free_nodal = read_xnodes(nodes.get_np(), basis, xnodes);
  }
  bool is_ok () const { return ok; }
  void eval (const Real& x, Real* const v) override {
    ::eval(nodes, get_xnodes(), x, v);
  }
  const Real* get_xnodes () const override {
    return free_nodal ? xnodes : islet::get_x_gll_special(nodes.get_np());
  }
  Int get_np () const override { return nodes.get_np(); }
};

} // namespace anon

extern "C" {

void offset_nodal_calc_xnodes_metrics (int np, int n, const int* subnp, const int* offst,
                                       double* metrics) {
  const auto nodes = make_offset_nodal(np, n, subnp, offst);
  calc_xnodes_metrics(nodes, islet::get_x_gll(np), metrics);
}

void calc_xnodes_metrics_from_basis_string (const char* basis, double* metrics) {
  Nodes nodes;
  const auto ok = nodes.init(basis);
  if ( ! ok) {
    printf("Invalid basis string: %s\n", basis);
    return;
  }
  Real xnodes[islet::np_max];
  const auto free_nodal = read_xnodes(nodes.get_np(), basis, xnodes);
  calc_xnodes_metrics(nodes,
                      free_nodal ? xnodes : islet::get_x_gll_special(nodes.get_np()),
                      metrics);
}

void calc_lebesgue_consts_from_basis_string (const char* basis, double* metrics) {
  Nodes nodes;
  const auto ok = nodes.init(basis);
  if ( ! ok) {
    printf("Invalid basis string: %s\n", basis);
    return;
  }
  Real xnodes[islet::np_max];
  const auto free_nodal = read_xnodes(nodes.get_np(), basis, xnodes);
  calc_lebesgue_consts(nodes,
                       free_nodal ? xnodes : islet::get_x_gll_special(nodes.get_np()),
                       metrics);
}

void run_thorough_diagnostics_from_basis_string (const char* basis) {
  {
    Real m[3];
    calc_xnodes_metrics_from_basis_string(basis, m);
    printf("    npm %1.4e %1.4e %1.4e\n", m[0], m[1], m[2]);
  }
  static const int ne_max = 11111;
  pum::Options po;
  po.threaded = true;
  po.ntrial = 33;
  po.mec_ne = 333;
  po.perturb = 0.01;
  printf("    ne,ndx_max %d po.ntrial %d po.mec_ne %d po.perturb %1.4f\n",
         ne_max, po.ntrial, po.mec_ne, po.perturb);
  auto b = std::make_shared<Basis>(basis);
  if ( ! b->is_ok()) return;
  {
    MaxEigComputer mec;
    const auto meam1 = mec.run(b->get_np(), ne_max, ne_max, 1e-13, true, b);
    printf("meam1 %1.4e\n", meam1);
  }
  {
    Real pum_max = 0;
    printf("pum:"); fflush(stdout);
    for (Int ne = 3; ne <= 15; ++ne) {
      po.ne = ne;
      pum::PerturbedUniformMeshMetric pum(b, po);
      const auto pum_val = pum.run();
      printf(" %1.1e", pum_val); fflush(stdout);
      pum_max = std::max(pum_max, pum_val);
    }
    printf("\npum_max %1.4e\n", pum_max);
  }
}

} // extern "C"
