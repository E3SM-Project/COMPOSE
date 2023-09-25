#include "islet_isl.hpp"
#include "islet_pum.hpp"

#include <cstdio>

static void write_array (const std::string& name, const std::vector<Real>& a) {
  printf("%s [", name.c_str());
  for (const auto& e : a) printf(" %1.5e", e);
  printf("]\n");
}

static void gen_plot_data (const islet::Operator::ConstPtr& op, const Int np,
                           const bool print_perturb = true) {
  static const Int nperturb = 48;
  std::vector<Real> perturb(nperturb), meam1(nperturb);
  for (Int i = 0; i < nperturb; ++i) perturb[i] = 0.1/std::pow(1.18, i);
  if (print_perturb) write_array(">> perturb", perturb);
  const auto oim = std::make_shared<islet::OperatorInterpMethod>(np, op);
  const auto im = std::make_shared<InterpMethod>(oim);
# pragma omp parallel for
  for (Int i = 0; i < nperturb; ++i) {
    meam1[i] = 0;
    pum::Options o;
    o.threaded = false;
    o.ntrial = 33;    
    o.perturb = perturb[i];
    for (const Int ne : {4, 7, 15})
      for (const Int mec_ne : {33}) {
        o.ne = ne;
        o.mec_ne = mec_ne;
        pum::PerturbedUniformMeshMetric pum(im, o);
        meam1[i] = std::max(meam1[i], pum.run());
      }
  }
  write_array(">> meam1", meam1);
}

int main (int argc, char** argv) {
  const auto gll_best = islet::Operator::create(islet::Operator::gll_best);
  const auto uofs = islet::Operator::create(islet::Operator::uniform_offset_nodal_subset);
  bool first = true;
  for (int np = 2; np <= islet::np_max; ++np) {
    if (np >= 4) {
      printf(">>> gll_best %d\n", np);
      gen_plot_data(gll_best, np, false);
    }
    printf(">>> uniform_offset_nodal_subset %d\n", np);
    gen_plot_data(uofs, np, first);
    first = false;
  }
}
