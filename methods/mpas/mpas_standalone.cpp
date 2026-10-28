#include "mpas_mpi.hpp"
#include "mpas_util.hpp"
#include "mpas_comm_util.hpp"
#include "mpas_partition.hpp"
#include "mpas_distribution.hpp"
#include "mpas_mesh_bulk_data.hpp"
#include "mpas_advect_ir_test.hpp"
// for unit tests
#include "mpas_siqp_test.hpp"
#include "mpas_lagrangian_levels.hpp"
#include "mpas_ci_util.hpp"
#include "mpas_accum.hpp"
#include "mpas_gallery.hpp"
#include "mpas_netcdf.hpp"

#include <iostream>

#ifdef MPAS_FPE
# include <xmmintrin.h>
#endif

static ir::OutputMgr::ConstPtr
setup_dev_output (const MpasMesh::ConstPtr& mm, const Comm::Ptr& he,
                  const ir::Params& p, const std::string& fprefix,
                  // Write initial and final points plus these.
                  const int n_internal_snaps = 11, const bool save_rho = false) {
  const auto om = std::make_shared<ir::OutputMgr>();
  om->time_step_write_stride = (p.nstep + n_internal_snaps)/(n_internal_snaps + 1);
  om->ncw = NetcdfWriter::create(he, mm, fprefix);
  if (p.vf->is_3d()) {
    if (p.vf->name() == "dcmip2012t11") {
      const Real zeta = 4900;
      const Int nlvl = mm->l.nmax;
      Int lvl;
      for (lvl = 0; lvl < nlvl; ++lvl)
        if (p.vf->ilvl_to_z(nlvl+1, lvl  ) <= zeta &&
            p.vf->ilvl_to_z(nlvl+1, lvl+1) >= zeta)
          break;
      const Real a = std::cos(5*(M_PI/6));
      const Vec axis = {std::sqrt(1 - square(a)), a, 0};
      for (int qi = 0; qi < p.ntracer; ++qi) {
        om->q_descs.push_back(om->ncw->add_cc_lonlat_time_slice("q", lvl, qi));
        om->q_descs.push_back(om->ncw->add_cc_vert_greatcircle_time_slice("q", axis, qi));
      }
    } else {
      const Int lvl = mm->l.nmax/2;
      om->q_descs.push_back(om->ncw->add_cc_lonlat_time_slice("q", lvl, 0));
      if (save_rho)
        om->rho_descs.push_back(om->ncw->add_cc_lonlat_time_slice("rho", lvl));
      if (p.vf->name() == "hadley") {
        const Vec axis = {1,1,0};
        om->q_descs.push_back(om->ncw->add_cc_vert_greatcircle_time_slice("q", axis, 0));
      } else {
        //om->q_descs.push_back(om->ncw->add_cc_vertlon_time_slice("q", 0., 0));
        const Vec axis = {0,0,1};
        om->q_descs.push_back(om->ncw->add_cc_vert_greatcircle_time_slice("q", axis, 0));
        if (save_rho)
          om->rho_descs.push_back(om->ncw->add_cc_vert_greatcircle_time_slice("rho", axis));
      }
    }
  } else {
    for (int qi = 0; qi < p.ntracer; ++qi)
      om->q_descs.push_back(om->ncw->add_cc_lonlat_time_slice("q", 0, qi));
  }
  om->ncw->end_definition();
  return om;
}

static bool cmp (const CellGraph& a, const CellGraph& b) {
  if (a.ncell() != b.ncell()) return false;
  const Int ncell = a.ncell();
  for (Int ci = 0; ci < ncell; ++ci)
    if (a.cells[ci] != b.cells[ci]) return false;
#if 1
  // The nbrs of a cell don't have to be ordered the same, but for now don't
  // permit that freedom.
  const Int nnbrs = a.cells[ncell];
  for (Int i = 0; i < nnbrs; ++i)
    if (a.nbrs[i] != b.nbrs[i]) return false;
#else
  // Or do.
  for (Int ci = 0; ci < ncell; ++ci) {
    const Int nnbr = a.cells[ci+1] - a.cells[ci];
    for (Int i = 0; i < nnbr; ++i) {
      const Int nbr = a.nbrs[a.cells[ci]+i];
      bool fnd = false;
      for (Int j = 0; j < nnbr; ++j)
        if (b.nbrs[b.cells[ci]+j] == nbr) {
          fnd = true;
          break;
        }
      if ( ! fnd) return false;
    }
  }
#endif
  return true;
}

struct Input {
  struct Command {
    enum Enum { run_ir, dev_io, dev_search, none };
    static Enum convert (const std::string& s) {
      if (s == "run_ir") return run_ir;
      if (s == "dev_io") return dev_io;
      if (s == "dev_search") return dev_search;
      return none;
    }
  };

  Command::Enum cmd;
  std::string part_fname, graph_fname, mesh_fname;
  std::string ir_opt_str;
  int nlevel, nstep, nsubstep, comm_nhalo, traj_nhalo;
  bool unittest, valid_input, ncout, lag_lvls_cell_int, zeta_eul_varying,
    perturb_dep_pts, bfb;

  Input(int argc, char** argv);
};

Input::Input (int argc, char** argv) {
  ir_opt_str = "na2elc";
  nlevel = 1;
  nstep = -1;
  comm_nhalo = 15;
  nsubstep = -1;
  traj_nhalo = -1;
  unittest = false;
  valid_input = true;
  ncout = false;
  lag_lvls_cell_int = true;
  zeta_eul_varying = false;
  perturb_dep_pts = false;
  bfb = true;
  int i = 1;
  if (argc > 1)
    cmd = Command::convert(argv[1]);
  if (cmd == Command::none)
    cmd = Command::run_ir;
  else
    ++i;
  for ( ; i < argc; ++i) {
    const std::string token(argv[i]);
    if (eq(token, "-p", "--part")) part_fname = argv[++i];
    else if (eq(token, "-m" , "--mesh")) mesh_fname = argv[++i];
    else if (eq(token, "-g" , "--graph")) graph_fname = argv[++i];
    else if (eq(token, "-o" , "--opt")) ir_opt_str = argv[++i];
    else if (eq(token, "-nlevel")) nlevel = std::atoi(argv[++i]);
    else if (eq(token, "-nstep")) nstep = std::atoi(argv[++i]);
    else if (eq(token, "-nsubstep")) nsubstep = std::atoi(argv[++i]);
    else if (eq(token, "-chalo")) comm_nhalo = std::atoi(argv[++i]);
    else if (eq(token, "-thalo")) traj_nhalo = std::atoi(argv[++i]);
    else if (eq(token, "-u", "--unittest")) unittest = true;
    else if (eq(token, "-ncout" )) ncout = true;
    else if (eq(token, "-llinterp")) lag_lvls_cell_int = false;
    else if (eq(token, "-zeulvar")) zeta_eul_varying = true;
    else if (eq(token, "-perturb")) perturb_dep_pts = true;
    else if (eq(token, "-bfboff")) bfb = false;
    else {
      std::cerr << "ERROR: Not an argument: " << token << "\n";
      valid_input = false;
    }
  }
}

static void print (const Input& in) {
  printf("mesh_fname = %s\n", in.mesh_fname.c_str());
  printf("graph_fname = %s\n", in.graph_fname.c_str());
  printf("part_fname = %s\n", in.part_fname.c_str());
  printf("ir_opt_str = %s\n", in.ir_opt_str.c_str());
  printf("nlevel = %d\n", in.nlevel);
  printf("nstep = %d\n", in.nstep);
  printf("comm_nhalo = %d\n", in.comm_nhalo);
  printf("traj_nhalo = %d\n", in.traj_nhalo);
  printf("unittest = %d\n", (int) in.unittest);
}

static void print_min_max_dzeta (const MpasMesh& mm, Comm& he) {
  if (mm.l.nmax < 2) return;
  Real lmin = mm.l.get_interface_z(0)[1] - mm.l.get_interface_z(0)[0],
    lmax = lmin, gmin, gmax;
  for (int ci = 0; ci < mm.t.ncell; ++ci) {
    const Real* zint = mm.l.get_interface_z(ci);
    for (int li = 0; li < mm.l.nmax; ++li) {
      const Real dz = zint[li+1] - zint[li];
      lmin = std::min(lmin, dz);
      lmax = std::max(lmax, dz);
    }
  }
  he.all_reduce(Comm::MIN, 1, &lmin, &gmin);
  he.all_reduce(Comm::MAX, 1, &lmax, &gmax);
  const Real q = gmin > 0 ? gmax/gmin : -1;
  if (he.amroot())
    printf("dzeta min %21.15e, max %21.15e, max/min %5.2f\n", gmin, gmax, q);
}

static void init_fpe () {
#ifdef MPAS_FPE
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() &
                         ~( _MM_MASK_INVALID |
                            _MM_MASK_DIV_ZERO |
                            _MM_MASK_OVERFLOW ));
#endif
}

// Example:
// $ gpmetis x1.2562.graph.info 3
// $ (np=3; f=x1.2562;
//    mpirun -np $np ./mpas_standalone -m $f.grid.nc -g $f.graph.info -p $f.graph.info.part.$np)
// -g is optional if -u (unit testing) is not requested.
int main (int argc, char** argv) {
  MPI_Init(&argc, &argv);
  const auto p = mpi::make_parallel(MPI_COMM_WORLD);

  Input in(argc, argv);
  if (p->amroot()) print(in);
  if ( ! in.ncout) init_fpe();

  int status = 0;
  do {
    if ( ! in.valid_input) break;
    const auto cg = create_CellGraph(in.mesh_fname);
    if ( ! cg) {
      if (p->amroot())
        printf("Could not read mesh_fname %s.\n", in.mesh_fname.c_str());
      break;
    }
    if (p->amroot()) printf("nbr cells %ld\n", cg->cells.size());
    const auto c2r = (p->size() > 1 ?
                      read_partition(in.part_fname) :
                      create_serial_partition(cg->ncell()));
    Int rmax = 0; for (const auto& r : c2r) rmax = std::max(rmax, r);
    if (rmax+1 != p->size()) {
      if (p->amroot())
        printf("Error: %d ranks but graph.info.part.%d\n", p->size(), rmax+1);
      status = -1;
      break;
    }

    const auto d = Distribution::create(p, c2r, *cg, in.comm_nhalo);
    assert(d);
    const auto mbd = MpasMeshBulkData
      ::create_and_update_distribution(in.mesh_fname, *d, in.nlevel);
    assert(mbd);
    const auto mm = mbd->get_mpas_mesh();
    const auto m = create_siqp_mesh(*mm);
    init_more_extras(*mm, *m);

    const auto he = Comm::create(d, mm);
    assert(he);
    if ( ! he->test()) printf("Comm::test failed\n");

    if (in.unittest) {
      if (p->amroot()) std::cout << "begin unit tests\n";
      if ( ! in.graph_fname.empty()) {
        const auto cg1 = read_cell_graph(in.graph_fname);
        if ( ! cmp(*cg1, *cg)) std::cout << "create_CellGraph FAIL\n";
      }
      if ( ! KahanBabushkaNeumaierSummer::unittest())
        std::cout << "KahanBabushkaNeumaierSummer::unittest FAIL\n";
      if ( ! order_invariant_accum::unittest_Accumulator())
        std::cout << "order_invariant_accum::unittest_Accumulator FAIL\n";
      if ( ! unittest_gallery(mm, *m, he))
        if (p->amroot()) std::cout << "unittest_gallery FAIL\n";
      if ( ! unittest_CellShapeFunction(mm, *m, he))
        if (p->amroot()) std::cout << "unittest_CellShapeFunction FAIL\n";
      if ( ! LagrangianLevels::unittest(mm, m, he, true))
        if (p->amroot()) std::cout << "LagrangianLevels::unittest FAIL\n";
      if ( Trajectory::test_search(mm, he) != 0)
        std::cout << "Trajectory::test_search FAIL\n";
      if (p->size() == 1) {
        if ( ! ir::unittest(mm, m, he, nullptr))
          std::cout << "ir::unittest FAIL\n";
        if ( ! test_siqp(he, mm)) std::cout << "test_siqp FAIL\n";
      }
      if (p->amroot()) std::cout << "end unit tests\n";
    }

    switch (in.cmd) {
    case Input::Command::run_ir: {
      ir::TestProblemOptions o;
      o.flow_ic_str = in.ir_opt_str;
      o.nstep = in.nstep;
      o.nsubstep = in.nsubstep;
      o.traj_nhalo = in.traj_nhalo;
      
      ir::Params p;
      p.lag_lvls_cell_int = in.lag_lvls_cell_int;
      p.perturb_dep_pts = in.perturb_dep_pts;
      p.bfb = in.bfb;
      ir::set(*mm, *he, o, p);
      
      if (in.zeta_eul_varying && p.vf->is_3d())
        mm->l.set_interface_spatially_varying_zeta_eul(
          mm->c, p.vf->ilvl_to_z(2, 0), p.vf->ilvl_to_z(2, 1));
      else
        set_interface_z_all_cols(*p.vf, *mm);
      init_zeta_rep(mm, *m, he);
      print_min_max_dzeta(*mm, *he);
      
      const auto op = ir::RemapFormOp::create(mm, m, he, p);
      ir::OutputMgr::ConstPtr om;
      if (in.ncout) om = setup_dev_output(mm, he, p, "tmp/ncwtest");
      ir::run_test_problem(op, o, om);
    } break;
    case Input::Command::dev_io: {
      NetcdfWriter::unittest(he, mm);
    } break;
    case Input::Command::dev_search: {
      Trajectory::test_search(mm, he, he->amroot());
    } break;
    case Input::Command::none: break;
    }
  } while (0);

  MPI_Finalize();
  return status;
}
