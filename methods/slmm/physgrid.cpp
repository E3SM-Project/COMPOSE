/* Study GLL <-> FV remap to support physgrid in E3SM.

   Options:
   -ne, -np, -nphys: Mesh params. (5, 4, 3)
   -nremap: Number of times to remap GLL -> FV -> GLL. (1)

     In this study, we assume both rho and q have to be remapped each
   direction. In practice, rho is not remapped from FV to GLL. So in that
   direction, ignore the remap operations for rho.

   Example:
       ./physgrid -ne 10 -np 4 -nphys 2 -nremap 1 -rho constant -q toychem1 -vis tmp/test
       (cd meas; hy physgrid.hy vis ../tmp/test ../fig/test)
 */

#include "slmm_gll.hpp"
#include "slmm_mesh.hpp"
#include "slmm_vis.hpp"

#include "slmmir_util.hpp"
#include "slmmir_physgrid.hpp"
using namespace pg;

struct Input {
  Int ne, np, nphys, nremap, res;
  Basis::Ptr basis;
  Limiter::Enum limiter;
  gallery::InitialCondition::Shape rho_ic, q_ic;
  bool print, omit_rho;
  std::string vis_fname_prefix;
  Fv2Gll::Type::Enum fv2gll_type;

  Input (int argc, char** argv)
    : ne(5), np(4), nphys(3), nremap(1), res(128), basis(std::make_shared<GLL>()),
      limiter(Limiter::caas), rho_ic(gallery::InitialCondition::XYZTrig),
      q_ic(gallery::InitialCondition::GaussianHills), print(false), omit_rho(false),
      fv2gll_type(Fv2Gll::Type::idem)
  {
    for (int i = 1; i < argc; ++i) {
      const std::string& token = argv[i];
      if (eq(token, "-ne"))
        ne = atoi(argv[++i]);
      else if (eq(token, "-np"))
        np = atoi(argv[++i]);
      else if (eq(token, "-nphys"))
        nphys = atoi(argv[++i]);
      else if (eq(token, "-nremap"))
        nremap = atoi(argv[++i]);
      else if (eq(token, "-res"))
        res = atoi(argv[++i]);
      else if (eq(token, "-rho"))
        rho_ic = gallery::InitialCondition::from_string(argv[++i]);
      else if (eq(token, "-q"))
        q_ic = gallery::InitialCondition::from_string(argv[++i]);
      else if (eq(token, "-lim"))
        limiter = Limiter::convert(argv[++i]);
      else if (eq(token, "-v", "-print"))
        print = true;
      else if (eq(token, "-omitrho"))
        omit_rho = true;
      else if (eq(token, "-vis"))
        vis_fname_prefix = argv[++i];
      else if (eq(token, "-basis"))
        basis = Basis::create(Basis::Type::convert(argv[++i]));
      else if (eq(token, "-fv2gll"))
        fv2gll_type = Fv2Gll::Type::convert(argv[++i]);
    }
  } 
};

struct Error {
  typedef std::shared_ptr<Error> Ptr;
  // GLL vs GLL
  Real l1_err, l2_err, linf_err, mass_cons, min_err, max_err;
  // GLL vs pg
  Real mass_cons_pg;
};

struct ProblemData {
  Array<Real> rho, q;

  void init(const Int n);

  void init(const Mesh& m, gallery::InitialCondition::Shape rho_ic,
            gallery::InitialCondition::Shape q_ic);
};

static Mesh::Ptr make_mesh (const Int ne, const Int np, const Basis::Ptr& basis) {
  Mesh::Ptr pm = std::make_shared<Mesh>();
  auto& m = *pm;
  m.np = np;
  m.basis = basis;
  m.tq_order = -1;
  mesh::make_cubedsphere_mesh(m.geo_p, m.geo_c2n, ne);
  mesh::make_cgll_from_geo(m.geo_p, m.geo_c2n, m.np, *m.basis, m.cgll_p, m.cgll_c2n);
  mesh::make_dgll_from_cgll(m.cgll_p, m.cgll_c2n, m.dglln2cglln, m.dgll_c2n);
  mesh::make_io_cgll_from_internal_cgll(m.cgll_p, m.cgll_c2n, m.cgll_io_c2n);
  fill_normals(m.geo_p, m.geo_c2n, m.geo_nml, m.geo_c2nml);
  return pm;
}

Error::Ptr measure_error (const PhysgridOps& ops, const MeshData& md,
                          const ProblemData& pd_gll0, const ProblemData& pd_pg,
                          const ProblemData& pd_gll1) {
  assert(static_cast<Int>(pd_gll0.q.size()) == md.ndgll);
  const auto err = std::make_shared<Error>();
  Real l1_num = 0, l1_den = 0, l2_num = 0, l2_den = 0, linf_num = 0, linf_den = 0;
  Real mass0 = 0, mass1 = 0;
  Real min0 = 1e3, max0 = -1e3, min1 = min0, max1 = max0;
  for (Int i = 0; i < md.ndgll; ++i) {
    const Real a = pd_gll0.q[i], b = pd_gll1.q[i];
    l1_num += std::abs(a - b)*md.spheremp[i];
    l1_den += std::abs(a)*md.spheremp[i];
    l2_num += square(a - b)*md.spheremp[i];
    l2_den += square(a)*md.spheremp[i];
    linf_num = std::max(linf_num, std::abs(a - b));
    linf_den = std::max(linf_den, std::abs(a));
    mass0 += pd_gll0.rho[i]*a*md.spheremp[i];
    mass1 += pd_gll1.rho[i]*b*md.spheremp[i];
    min0 = std::min(min0, a); max0 = std::max(max0, a);
    min1 = std::min(min1, b); max1 = std::max(max1, b);
  }
  err->l1_err = l1_num/l1_den;
  err->l2_err = std::sqrt(l2_num/l2_den);
  err->linf_err = linf_num/linf_den;
  err->mass_cons = (mass1 - mass0)/mass0;
  err->min_err = std::min(min1 - min0, 0.0);
  err->max_err = std::max(max1 - max0, 0.0);
  Real mass_pg = 0;
  const Real wt_pg = 1.0/square(ops.nphys);
  for (Int i = 0, n = md.fv_metdet.size(); i < n; ++i)
    mass_pg += md.fv_metdet[i]*wt_pg*pd_pg.rho[i]*pd_pg.q[i];
  err->mass_cons_pg = (mass_pg - mass0)/mass0;
  return err;
}

void print (const Input& in, const Error& e) {
  printf("dpg> rho %s q %s\n",
         gallery::InitialCondition::to_string(in.rho_ic),
         gallery::InitialCondition::to_string(in.q_ic));
  printf("dpg> ne %d np %d nphys %d nremap %d fv2gll %s lim %s\n",
         in.ne, in.np, in.nphys, in.nremap,
         Fv2Gll::Type::convert(in.fv2gll_type).c_str(),
         Limiter::convert(in.limiter).c_str());
  const auto eps = std::numeric_limits<Real>::epsilon();
  const std::string msg =
    ((in.limiter != Limiter::none && (e.min_err < -eps || e.max_err > eps)) ||
     std::abs(e.mass_cons) > 1e3*eps) ?
    " ERROR" : "";
  printf("dpg> l1 %8.2e l2 %8.2e linf %8.2e cons %9.2e %9.2e min %9.2e max %8.2e%s\n",
         e.l1_err, e.l2_err, e.linf_err, e.mass_cons, e.mass_cons_pg, e.min_err,
         e.max_err, msg.c_str());
}

void ProblemData::init (const Mesh& m, gallery::InitialCondition::Shape rho_ic,
                        gallery::InitialCondition::Shape q_ic) {
  const auto dnn = nslices(m.dglln2cglln);
  init(dnn);
  for (Int i = 0; i < dnn; ++i) {
    const auto k = m.dglln2cglln[i];
    const auto n = slice(m.cgll_p, k);
    Real lat, lon;
    xyz2ll(n[0], n[1], n[2], lat, lon);
    gallery::InitialCondition::init(rho_ic, 1, &lat, &lon, &rho[i]);
    // Keep rho well behaved, as some ICs go to 0. We're not interested in the
    // details of positivity of rho in this study.
    rho[i] += 1;
    gallery::InitialCondition::init(q_ic, 1, &lat, &lon, &q[i]);
  }
}

void ProblemData::init (const Int n) {
  rho.optclear_and_resize(n, 0);
  q.optclear_and_resize(n, 0);
}

void remap (const Mesh& m, const MeshData& md, const Gll2Fv& pg,
            const ProblemData& d_dyn, ProblemData& d_phys,
            const Limiter::Enum limiter) {
  const Int nphys = pg.get_nphys(), nphys2 = nphys*nphys,
    np = pg.get_np(), np2 = np*np;
  assert(static_cast<Int>(d_dyn.q.size()) == md.nelem * np2);
  assert(static_cast<Int>(d_phys.q.size()) == md.nelem * nphys2);
  for (Int ie = 0; ie < md.nelem; ++ie)
    pg.remap(&md.gll_metdet[ie*np2], &md.fv_metdet[ie*nphys2], limiter,
             &d_dyn.rho[ie*np2], &d_dyn.q[ie*np2],
             &d_phys.rho[ie*nphys2], &d_phys.q[ie*nphys2]);
}

void remap (const Mesh& m, const MeshData& md, const Fv2Gll& pg,
            const ProblemData& d_phys, ProblemData& d_dyn,
            const Limiter::Enum limiter, const bool omit_rho = false) {
  const Int nphys = pg.get_nphys(), nphys2 = nphys*nphys,
    np = pg.get_np(), np2 = np*np;
  const bool limit = nphys > 1 && limiter != Limiter::none;
  assert(static_cast<Int>(d_dyn.q.size()) == md.nelem * np2);
  assert(static_cast<Int>(d_phys.q.size()) == md.nelem * nphys2);
  Array<Real> qlos, qhis;
  if (limit) {
    // Obtain limiter bounds. This will use Homme's limiter halo exchange.
    qlos.optclear_and_resize(md.nelem);
    qhis.optclear_and_resize(md.nelem);
    for (Int ie = 0; ie < md.nelem; ++ie) {
      const Real* q = &d_phys.q[ie*nphys2];
      Real qlo = q[0], qhi = q[0];
      for (Int i = 1; i < nphys2; ++i) qlo = std::min(qlo, q[i]);
      for (Int i = 1; i < nphys2; ++i) qhi = std::max(qhi, q[i]);
      qlos[ie] = qlo;
      qhis[ie] = qhi;
    }
  }
  // Remap rho, q = Q/rho.
  for (Int ie = 0; ie < md.nelem; ++ie) {
    Real qlo = -1, qhi = -1;
    if (limit) {
      // Bounds from cell ie and its neighbors.
      qlo = qlos[ie];
      qhi = qhis[ie];
      for (Int j = md.geo_c2cnbrs_ptr[ie]; j < md.geo_c2cnbrs_ptr[ie+1]; ++j) {
        qlo = std::min(qlo, qlos[md.geo_c2cnbrs[j]]);
        qhi = std::max(qhi, qhis[md.geo_c2cnbrs[j]]);
      }
    }
    pg.remap(&md.gll_metdet[ie*np2], &md.fv_metdet[ie*nphys2],
             limiter, qlo, qhi,
             &d_phys.rho[ie*nphys2], &d_phys.q[ie*nphys2],
             &d_dyn.rho[ie*np2], &d_dyn.q[ie*np2], ! omit_rho);
  }
  if ( ! omit_rho) {
    // q -> Q
    for (Int i = 0; i < md.ndgll; ++i)
      d_dyn.q[i] *= d_dyn.rho[i];
    // DSS rho, Q.
    md.d2cer->dss(d_dyn.rho.data(), md.wrk.data());
    md.d2cer->dss(d_dyn.q.data(), md.wrk.data());
    // Q -> q
    for (Int i = 0; i < md.ndgll; ++i)
      d_dyn.q[i] /= d_dyn.rho[i];
  } else {
    md.d2cer->dss_q(d_dyn.rho.data(), d_dyn.q.data(), md.wrk.data());
  }
  if (nphys == 1)
    for (Int ie = 0; ie < md.nelem; ++ie)
      pg.reconstruct_nphys1(&md.gll_metdet[ie*np2],
                            &d_dyn.rho[ie*np2], &d_dyn.q[ie*np2]);
}

int mods = 0, calls = 0;

static void run (const Input& in) {
  const auto mesh = make_mesh(in.ne, in.np, in.basis);
  const auto ops = std::make_shared<PhysgridOps>(*mesh, in.np, in.nphys, in.fv2gll_type);
  const auto pdgll0 = std::make_shared<ProblemData>();
  pdgll0->init(*mesh, in.rho_ic, in.q_ic);
  const auto pdphys = std::make_shared<ProblemData>();
  pdphys->init(ops->mesh_data->nelem * in.nphys * in.nphys);
  const auto pdgll = std::make_shared<ProblemData>(*pdgll0);
  //pdgll->init(pdgll0->rho.size());
  if (in.print) {
    ops->gll2fv->print();
    ops->fv2gll->print();
  }
  vis::VisWriter::Ptr iw_gll, iw_pg;
  Array<Real> wrk;
  if ( ! in.vis_fname_prefix.empty()) {
    const Int res = in.res;
    const auto& m = *mesh;
    const auto op_gll = std::make_shared<vis::BilinGLLToLatLon>(
      m.cgll_p, m.cgll_io_c2n, 2*res+1, 4*res+1);
    iw_gll = std::make_shared<vis::VisWriter>(op_gll, in.vis_fname_prefix + "_gll.bin");
    iw_pg = make_pg_vis(in.vis_fname_prefix, in.ne, in.nphys, res);
    wrk.optclear_and_resize(pdgll0->q.size());
  }
  ProblemData pdgllc(*pdgll0);
  ProblemData* pdglls[] = {pdgll.get(), &pdgllc};
  const auto& md = *ops->mesh_data;
  for (Int it = 0; it < in.nremap; ++it) {
    std::swap(pdglls[0], pdglls[1]);
    //pr("d->p");
    if (it == 0 && iw_gll) {
      if ( ! in.omit_rho) {
        md.d2cer->d2c(pdglls[0]->rho.data(), wrk.data(), true);
        iw_gll->write(wrk.data());
      }
      md.d2cer->d2c(pdglls[0]->q.data(), wrk.data(), true);
      iw_gll->write(wrk.data());
    }
    remap(*mesh, md, *ops->gll2fv, *pdglls[0], *pdphys, in.limiter);
    if (iw_pg) {
      if ( ! in.omit_rho) iw_pg->write(pdphys->rho.data());
      iw_pg->write(pdphys->q.data());
    }
    //printf("%5d %5d\n", mods, calls); mods = calls = 0; pr("p->d");
    remap(*mesh, md, *ops->fv2gll, *pdphys, *pdglls[1], in.limiter, in.omit_rho);
    if (iw_gll) {
      if ( ! in.omit_rho) {
        md.d2cer->d2c(pdglls[1]->rho.data(), wrk.data(), true);
        iw_gll->write(wrk.data());
      }
      md.d2cer->d2c(pdglls[1]->q.data(), wrk.data(), true);
      iw_gll->write(wrk.data());
    }
    //printf("%5d %5d\n",mods,calls);
  }
  const auto err = measure_error(*ops, md, *pdgll0, *pdphys, *pdglls[1]);
  print(in, *err);
}

int main (int argc, char** argv) {
  Kokkos::initialize(argc, argv); {
    Input in(argc, argv);
    run(in);
  } Kokkos::finalize_all();
}
