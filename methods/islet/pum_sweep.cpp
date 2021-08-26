#include <cstdlib>

#include "islet_isl.hpp"
#include "islet_pum.hpp"

// Illustrate why a 1-element hop is the key thing to study. Report meam1 at x
// in [-1,1] if meam1 >= tol. x, meam1 are 1-1.
static void run (const islet::Operator::ConstPtr& op,
                 const Int np, const Int nx, const Int ntrial,
                 std::map<Real,Real>& dx2meam1) {
  const auto oim = std::make_shared<islet::OperatorInterpMethod>(np, op);
  printf("%s\n", op->get_basis_string(np).c_str());
  const auto x = op->get_xnodes(np);
  printf("x:");
  for (Int i = np/2; i < np; ++i)
    printf(" %7.5f", x[i]);
  printf("\n");
  pum::Options o;
  o.threaded = true;
  o.perturb = 0.01;
  for (const Int mec_ne : {3, 33, 333})
    for (const Int ne : {4, 7, 15}) {
      printf("ntrial %d mec_ne %d ne %d\n", ntrial, mec_ne, ne);
      o.ntrial = ntrial;
      o.mec_ne = mec_ne;
      o.ne = ne;
      pum::PerturbedUniformMeshMetric pum(oim, o);
      pum.sweep_and_collect_amplitudes(nx, 1e-13, dx2meam1, false);
    }
  printf("final\n");
  for (const auto& e : dx2meam1)
    printf("dx %1.16e meam1 %1.16e\n", e.first, e.second);
}

int main (int argc, char** argv) {
  if (argc < 5) {
    printf("%s np nx ntrial (0 - natural, 1 - gll_best, 2 - uniform)\n", argv[0]);
    return -1;
  }
  const Int np = std::atoi(argv[1]);
  const Int nx = std::atoi(argv[2]);
  const Int ntrial = std::atoi(argv[3]);
  const Int opcode = std::atoi(argv[4]);
  if (np < 4 || nx < 2 || ntrial < 1 || opcode < 0 || opcode > 2) {
    printf("bad input");
    return -1;
  }
  const auto op = islet::Operator::create(opcode == 0 ?
                                          islet::Operator::gll_natural :
                                          opcode == 1 ?
                                          islet::Operator::gll_best :
                                          islet::Operator::uniform_offset_nodal_subset);
  std::map<Real,Real> dx2meam1;
  run(op, np, nx, ntrial, dx2meam1);
  return 0;
}
