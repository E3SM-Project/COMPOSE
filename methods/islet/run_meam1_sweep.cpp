#include "islet_isl.hpp"
#include "islet_maxeigcomp.hpp"
#include "islet_pum.hpp"

#include <cstdio>

static void
run_sweep (const islet::Operator::ConstPtr& op, const Int np) {
  static const Int nx = 192;
  std::vector<Real> dxs(nx), meam1s(nx);
  const auto oim = std::make_shared<islet::OperatorInterpMethod>(np, op);
  const auto im = std::make_shared<InterpMethod>(oim);
# pragma omp parallel for schedule(static,1)
  for (Int ix = 0; ix < nx; ++ix) {
    MaxEigComputer mec(false);
    const Real dx = 0.5*(Real(ix+1)/nx);
    dxs[ix] = dx;
    Real mea;
    mec.compute(*im, dx, 1024, &mea);
    meam1s[ix] = mea - 1;
  }
  for (Int ix = 0; ix < nx; ++ix)
    printf("%23.15e %23.15e\n", dxs[ix], meam1s[ix]);
}

int main (int argc, char** argv) {
  const int np = argc == 2 ? std::atoi(argv[1]) : -1;
  if (argc != 2 || np == -1) {
    printf("%s np\n", argv[0]);
    return -1;
  }
  {
    const auto gll_natural = islet::Operator::create(islet::Operator::gll_natural);
    printf("gll_natural %d\n", np);
    run_sweep(gll_natural, np);
  }
  {
    const auto gll_best = islet::Operator::create(islet::Operator::gll_best);
    printf("gll_best %d\n", np);
    run_sweep(gll_best, np);
  }
  {
    const auto uofs = islet::Operator::create(islet::Operator::uniform_offset_nodal_subset);
    printf("uniform_offset_nodal_subset %d\n", np);
    run_sweep(uofs, np);
  }
  return 0;
}
