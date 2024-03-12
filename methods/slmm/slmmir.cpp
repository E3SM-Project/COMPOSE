/* SLMMIR: Semi-Lagrangian Multi-Moment Incremental-Remap
           algorithm analyzer and tester.

   This program runs our various SLMMIR algorithms on test cases. There are a
   number of options, as follows:

   -method: Integration method, one of {ir, cdg, isl, pisl}. Default: ir. ir is
     incremental remap; density Q_ki is propagated and then remapped. cdg is
     characteristic discontinuous Galerkin. csl is classical semi-Lagrangian,
     but with rho remapped; this is a technical option and in practice should
     not be used. pcsl is pure ISL and should be used. pcsl is super fast
     compared with the others, so use it if you're just wanting to get a
     solution of some sort. Alias pairs: (isl,csl), (pisl,pcsl).
   -ode: The ODE test case to run, one of {dcmip1d3ll, nondivergent, divergent,
     rotate, movingvortices}. Default: divergent.
   -ic: Initial condition, in {xyztrig, gaussianhills, cosinebells,
     slottedcylinders, correlatedcosinebells, vortextracer}. Default:
     xyztrig. This option can be specified multiple times to run several fields
     together.
   -xyz: Integrate in (x,y,z) space instead of lat-lon.
   -T: Days to integrate. Default: 12, which is the repeat period of the flow
     field.
   -nsteps: Number of time steps per 12 days. Default: 120.
   -timeint: Time integration method in {exact, line, interp, interpline}.
    Default: exact. (Exact is not truly exact, but the time integration is done
    sufficiently accurately as not to be part of the error.)  The 'interp'
    options integrate on the np4 mesh and interpolate to the np>4 tracer mesh.
   -mesh: Mesh type. Default: geometric. gllsubcell and runisubcell are
    experimental.
   -ne: Number of elements per side of the square. Default: 5.
   -nonunimesh: 0 (default) for quasiuniform.
   -np: Spectral element order, in {2, 3, 4}. Default: 4. For pcsl, np=4 uses a
     modified basis for stability.
   -pg: Use physgrid for physics. 0 means none.
   -tq: Triangle quadrature order, in {4, 8, 12, 14, 20}: Default: 20 for np=4,
     12 otherwise.
   -dmc: Discrete mass conservation method, in {equality-sphere (es),
     equality-homme (eh), facet (f), equality-facet (ef), none (n),
     global-equality-homme (geh)}. Default: none.
   -mono: Apply a shape preservation filter, in {none, qlt, qlt-pve, caas,
     caas-pve, mn2, caas-node}. This applies to the global problem. none - no
     limiter, qlt - tree-based limiter, caas - ClipAndAssuredSum, mn2 - 2-norm
     minimization (not threaded), *-pve: fix mass subject to >= 0. Default:
     none. caas-node runs CAAS on each node, so the element-local limiter is not
     used.
   -lim: mn2 (min-2-norm), caas, caags, none. Default: mn2. Only active if -mono
     is not none. This applies to the cell-local problem.
   -subcellbounds: If set, use only the 4 surrounding points for bounds.
   -o: A string that is a prefix for all output files. Default: tmp/out
   -we: Output write interval, in time steps. Default: 1. Set to 0 to turn off
     output.
   -io-type: I/O output type, in {netcdf, internal}. Default: netcdf. If netcdf,
     write output for use in paraview; if internal, latlon output for python.
   -io-nodss: In the I/O DSS for q, inject just one value rather than summing
     according to the DSS.
   -io-start-cycle: Wait to start output, except ICs, until this cycle. Cycle
    count starts at 1.
   -io-recon: For the internal I/O type, choose the reconstruction from {bilin,
    const}. 'bilin' and 'const' act on subcells. 'bilin' is bilinear interp;
    'const' is a simple average of nodal values. Default: bilin.
   -res: Output resolution; depends on writer. If internal, half number of
     latitude points. Default: 64.
   -rit: Record scalar measurements in time. Output is a matlab file with a
     struct that can be loaded. The file is somewhat self-documenting.
   -lauritzen: Run Lauritzen diagnostics. Only active if -rit.
   -lauritzen-io: Include expensive I/O.
   -footprint: Track ISL and CISL footprints.
   -midpoint-check: Output error analysis for midpoint of run w.r.t. to a
    high-res 1-step solution.
   -rotate-grid: If provided, rotate the grid so that, in particular, it's not
    aligned with the N-S poles and other features of the flows.
   -rhot0: If provided, freeze rho at its initial value. Impl'ed for only CSL
    methods.

   Debug and analysis options:

   --write-binary-at: Write a potentially huge binary file at a given time step,
     to be read by msik('mc_read_slmmir_binary').

   Some notes:

   -dmc eh uses QOS and applies a cell-local equality constraint to match Homme
   mass in the cell; this drops the OOA of the shape functions by 1. geh applies
   a single constraint to the whole system and does not drop the OOA, but it
   won't work with shape preservation. geh was done before QLT was created. With
   QLT, one can do the equivalent of geh plus do shape preservation. However,
   this combination is not implemented.

   -mono vs -lim. There are two parts to shape preservation, even w/o the
   problem of tracer consistency. First is that cell mass (0th moment) must be
   in bounds. For p>1, this is not assured from the discretization alone. If
   -mono is not none, then mass is redistributed among cells to get all cell
   masses in bounds. This is the purpose of QLT (-mono qlt).
     Once cell masses are in bounds, then any number of methods can be used
   within each cell to limit the higher moments. The default is -lim mn2, which
   is the method of Guba et al JCP 2014. Algorithms for the global problem can
   also be applied locally.
     In practice (tracer transport coupled to another method that solves the
   continuity equation among others), tracer consistency is also an issue. -mono
   takes care of this, too.

   -timeint = interp explores whether it makes sense to try a finer tracer grid
   relative to the dynamics (in this case, imposed velocity) grid. timeint_v_np
   could be made an option, but for now is hardcoded to the value of
   interest. It is 4 for np >= 4 and np for np <= 3.
     Pass -prefine n, (default 0) to indicate a particular p-refinement
   experiment impl.
     0 - Default. Run the problem as usual, but possibly with -timeint interp*.
   Subsequent options implicitly have -timeint interp*.
     1 - Difference from 0 is that rho is integrated on the v-grid.
   Subsequent options move all data back to the v-grid between time steps;
   toychem runs on the v-grid.
     5 - Interp between np=4 v-grid and p-refined t-grid. Diagnostics are on
   v-grid.
     -basis defaults to "GllNodal". Other options are "GllOffsetNodal",
   "UniformOffsetNodal" and a basis string.
 */

#include "slmm_defs.hpp"
#include "slmm_mesh.hpp"
#include "slmm_spf.hpp"
#include "slmm_gll.hpp"
#include "slmm_io.hpp"
#include "slmm_nla.hpp"
#include "slmm_gallery.hpp"
#include "slmm_accum.hpp"
#include "slmm_debug.hpp"
#include "slmm_fit_extremum.hpp"
#include "slmm_util.hpp"
#include "slmm_vis.hpp"
using namespace slmm;

#include "slmmir_util.hpp"
#include "slmmir_mesh.hpp"
#include "slmmir_remap_data.hpp"
#include "slmmir_mono_data.hpp"
#include "slmmir_remapper.hpp"
#include "slmmir_d2c.hpp"
#include "slmmir_time_int.hpp"
#include "slmmir_p_refine.hpp"
#include "slmmir_lauritzen_diag.hpp"
#include "slmmir_snapshot.hpp"

#include <fstream>

// Some debug and code stuff.
namespace {
class Debug {
  int index_;
  std::string filename_;
  bool on_;

public:
  Debug ()
    : index_(1), filename_("dbgout.m"), on_(true)
  {
#ifdef SLMM_DEBUG
    FILE* fid = fopen(filename_.c_str(), "w");
    fclose(fid);
#endif
  }

  void advance () { ++index_; }

  void set_on (const bool set) { on_ = set; }

  template <typename CV3s>
  void write_p (const std::string& name, const CV3s& p) {
#ifdef SLMM_DEBUG
    if ( ! on_) return;
    FILE* fid = fopen(filename_.c_str(), "a");
    fprintf(fid, "%s{%d} = [", name.c_str(), index_);
    for (Int ip = 0; ip < nslices(p); ++ip)
      fprintf(fid, " %1.15e %1.15e %1.15e;", p(ip,0), p(ip,1), p(ip,2));
    fprintf(fid, "].';\n");
    fclose(fid);
#endif
  }

  template <typename CIs>
  void write_c2n (const std::string& name, const CIs& e) {
#ifdef SLMM_DEBUG
    if ( ! on_) return;
    FILE* fid = fopen(filename_.c_str(), "a");
    fprintf(fid, "%s{%d} = [", name.c_str(), index_);
    for (Int ie = 0; ie < nslices(e); ++ie) {
      for (Int k = 0; k < szslice(e); ++k)
        fprintf(fid, " %d", e(ie,k)+1);
      fprintf(fid, ";");
    }
    fprintf(fid, "].';\n");
    fclose(fid);
#endif
  }

  void write (const std::string& name, const BlockMatrix<Real, Size, Int>& m) {
#ifdef SLMM_DEBUG
    if ( ! on_) return;
    FILE* fid = fopen(filename_.c_str(), "a");
    fprintf(fid, "tmp = [");
    const Size* rowptr = m.rowptr();
    const Int* colidx = m.colidx();
    for (Int R = 0; R < m.M(); ++R)
      for (Int J = rowptr[R]; J < rowptr[R+1]; ++J) {
        const Int C = colidx[J];
        const Real* const block = m.block(R, C);
        for (Int r = 0, k = 0; r < m.m(); ++r)
          for (Int c = 0; c < m.n(); ++c, ++k)
            fprintf(fid, "%d %d %1.15e\n", m.m()*R + r + 1,
                    m.n()*C + c + 1, block[k]);
      }
    fprintf(fid, "];\n");
    fprintf(fid, "%s{%d} = sparse(tmp(:,1),tmp(:,2),tmp(:,3),%d,%d);\n",
            name.c_str(), index_, m.M()*m.m(), m.N()*m.n());
    fclose(fid);
#endif
  }

#ifdef SLMM_DEBUG // not used
  void write (const std::string& name, const Real* const a, const Int n) {
    if ( ! on_) return;
    FILE* fid = fopen(filename_.c_str(), "a");
    fprintf(fid, "%s{%d} = [", name.c_str(), index_);
    for (Int i = 0; i < n; ++i)
      fprintf(fid, " %1.15e", a[i]);
    fprintf(fid, "].';\n");
    fclose(fid);
  }
#endif
};
static Debug gdbg;

static bool write (std::ofstream& os, const FullMassMatrix::MT& a) {
  using slmm::write;
  const auto nnz = a.rowptr()[a.M()];
  return (write(os, a.M()) &&
          write(os, a.N()) &&
          write(os, a.m()) &&
          write(os, a.n()) &&
          write(os, a.M() + 1, a.rowptr()) &&
          write(os, nnz, a.colidx()) &&
          write(os, nnz*a.m()*a.n(), a.blockrow(0)));
}

// For use by msik('mc_read_slmmir_binary').
static void write_binary (
  const FullMassMatrix::MT& M_fac, const RemapData::MT& T,
  const Real* const dgbfit, const Real* const dgbfit_gll,
  const Real* const Jt, const Real* const Js,
  const Int* const dglln2cglln, const std::string& out_fn)
{
  using slmm::write;
  std::ofstream os(out_fn.c_str(), std::ofstream::binary);
  if ( ! os.is_open())
    throw std::runtime_error("write_binary: could not open file.");
  const Int n = T.M()*T.m();
  Int nc = 0;
  for (Int i = 0; i < n; ++i) nc = std::max(nc, dglln2cglln[i]);
  ++nc;
  if ( ! (write(os, T) &&
          write(os, M_fac) &&
          write(os, n, dgbfit) &&
          write(os, n, dgbfit_gll) &&
          write(os, n, Jt) &&
          write(os, n, Js) &&
          write(os, n, dglln2cglln)))
    throw std::runtime_error("write_binary: could not write file.");
}
} // anon namespace

#ifdef SLMM_TIME
timeval Timer::t_start_[Timer::NTIMERS];
double Timer::et_[Timer::NTIMERS];
#endif

struct IntegrateOptions {
  // If fwd, integrate mesh forward in time from t_{n-1} to t_n.
  enum Direction { fwd, bwd };
  Direction stepping;
  // Incremental remap or classical semi-Lagrangian.
  Method::Enum method;
  // Each step, and in error measurement, convert dgll <-> cgll.
  bool d2c;
  // Apply a shape-preserving filter.
  Filter::Enum filter;
  // Cell-local limiter to use in filter.
  Limiter::Enum limiter;
  // Use only the 4 surrounding points for bounds.
  bool subcell_bounds;
  // Relax bounds using a quadratic fit to the cell field.
  bool fitext;
  // Record assessment quantities at each time step.
  bool record_in_time;
  // Assess error at day 6 against higher-res, 1-step reference.
  bool check_midpoint;
  // Run Lauritzen et al GMD 2012 diagnostics, if record_in_time.
  bool lauritzen_diag, lauritzen_diag_io;
  // Test tracer consistency by randomly perturbing rho.
  Real perturb_rho;
  // Discrete mass conservation.
  Dmc::Enum dmc;
  // Track ISL and CISL footprints.
  bool track_footprint;
  // Netcdf or custom.
  IOType::Enum io_type;
  vis::Reconstruction::Enum io_recon;
  bool io_no_dss;
  Int io_start_cycle;
  // If I/O is internal type, lat,lon grid resolution parameter.
  Int output_resolution;
  Int pg;
  bool rhot0, vortex_problem;
};

struct Input {
  struct InitialCondition {
    std::string name;
    Int n; // Repeat this IC n times (for performance profiling).
    InitialCondition (const std::string& name, const Int& n)
      : name(name), n(n) {}
  };
  std::string output_fn, ode, program_name;
  std::vector<InitialCondition> initial_conditions;
  Real T;
  TimeInt::Enum timeint;
  Int timeint_v_np, p_refine_experiment;
  Int ne, nonunimesh, nsteps, write_every, np, tq_order;
  MeshType::Enum mesh_type;
  bool xyz_form; // Integrate in (x,y,z) space instead of (lat,lon).
  IntegrateOptions integrate_options;
  std::string basis;
  bool rotate_grid, ode_tight;

  // Super-duper debug/analysis options.
  struct Debug {
    Int write_binary_at;
    Debug () : write_binary_at(-1) {}
  };
  Debug debug;

  Input(Int argc, char** argv);
  void print(std::ostream& os) const;
};

// _s is start and _e is end.
struct Output {
  Real
    l1_err, l2_err, max_err, // error norm of q
    mass_s, mass_e, // mass of Q at start, end
    min_s, max_s, min_e, max_e, // min, max of q at start, end
    et_timestep, mem_hightwater,
    mass_gll_s, mass_gll_e; // mass of Q using Homme def
};

static void print_error (
  const Mesh& m, const ARealArray& F_gll,
  const ARealArray& F_sphere,
  const Real* const fs, const Real* const ds,
  const Real* const fe, const Real* const de, Output& out)
{
  Real l2_num = 0, l2_den = 0, max_num = 0, max_den = 0;
  out.max_s = -1e300; out.min_s = 1e300;
  out.max_e = -1e300; out.min_e = 1e300;
  out.mass_s = 0; out.mass_e = 0;
  out.mass_gll_s = 0; out.mass_gll_e = 0;
  siqk::TriangleQuadrature tq;
  siqk::RawConstVec3s tq_bary;
  siqk::RawConstArray tq_w;
  tq.get_coef(m.tq_order, tq_bary, tq_w);
  const auto& c2n = m.geo_c2n;
  const Int np = m.np, np2 = square(np);
  // GLL mass conservation.
  for (Int ci = 0; ci < nslices(c2n); ++ci)
    for (Int j = 0, basis_idx = 0; j < m.np; ++j)
      for (Int i = 0; i < m.np; ++i, ++basis_idx) {
        const Int k = ci*np2 + basis_idx;
        const Real w = F_gll[k];
        const Int idx = k;
        out.mass_gll_s += w * fs[idx]*ds[idx];
        out.mass_gll_e += w * fe[idx]*de[idx];
      }
  // Mass conservation wrt quadrature approximation of exact integrals.
  for (Int ci = 0; ci < nslices(c2n); ++ci)
    for (Int j = 0, basis_idx = 0; j < m.np; ++j)
      for (Int i = 0; i < m.np; ++i, ++basis_idx) {
        const Int k = ci*np2 + basis_idx;
        const Real w = F_sphere[k];
        const Int idx = k;
        const Real e = fe[idx] - fs[idx];
        out.mass_s += w * fs[idx]*ds[idx];
        out.mass_e += w * fe[idx]*de[idx];
        l2_num += w * square(e);
        l2_den += w * square(fs[idx]);
        max_num = std::max(max_num, std::abs(e));
        max_den = std::max(max_den, std::abs(fs[idx]));
        out.min_s = std::min(out.min_s, fs[idx]);
        out.max_s = std::max(out.max_s, fs[idx]);
        out.min_e = std::min(out.min_e, fe[idx]);
        out.max_e = std::max(out.max_e, fe[idx]);
      }
  out.l2_err = std::sqrt(l2_num/l2_den);
  out.max_err = max_num/max_den;
  printf("> re l2 %9.3e max %9.3e\n", out.l2_err, out.max_err);
  // These are mass conservation measurements of the final field w.r.t. the
  // initial condition. In contrast, the "C" output is the max over the
  // per-time-step mass conservation error.
  printf("> [cv] re %10.3e\n", reldif(out.mass_s, out.mass_e));
  printf("> [cv gll] re %10.3e\n", reldif(out.mass_gll_s, out.mass_gll_e));
  printf("> [mo] min %10.3e %10.3e [%10.3e] max %10.3e %10.3e [%10.3e]\n",
         out.min_s, out.min_e, out.min_e - out.min_s,
         out.max_s, out.max_e, out.max_e - out.max_s);
}

static void print_one_liner (const Input& in, const Output& out) {
  std::cout << "<OL> method " << in.integrate_options.stepping
            << " ode " << in.ode << " ic " << in.initial_conditions[0].name
            << " T " << in.T << " np " << in.np << " ne " << in.ne
            << " nonunimesh " << in.nonunimesh << " tq " << in.tq_order
            << " nsteps " << in.nsteps << " mono " << in.integrate_options.filter;
  printf(" re l2 %9.3e max %9.3e", out.l2_err, out.max_err);
  printf(" cv re %9.3e", reldif(out.mass_s, out.mass_e));
  printf(" cvgll re %9.3e", reldif(out.mass_gll_s, out.mass_gll_e));
  printf(" mo min %9.3e %9.3e %9.3e max %9.3e %9.3e %9.3e",
         out.min_s, out.min_e, out.min_e - out.min_s,
         out.max_s, out.max_e, out.max_e - out.max_s);
  printf(" et ts %9.3e nthr %d", out.et_timestep, omp_get_max_threads());
  std::cout << " prog " << in.program_name;
  std::cout << " xyz " << in.xyz_form;
  std::cout << " d2c " << in.integrate_options.d2c;
  std::cout << " dmc " << Dmc::convert(in.integrate_options.dmc);
  std::cout << " method " << Method::convert(in.integrate_options.method);
  std::cout << " limiter " << Limiter::convert(in.integrate_options.limiter);
  if (in.integrate_options.subcell_bounds)
    std::cout << " subcellbounds";
#ifdef RELAX_TIME
  std::cout << " RELAX_TIME";
#endif
  std::cout << " timeint " << TimeInt::convert(in.timeint) << " " << in.timeint_v_np;
  std::cout << "\n";
}

static void init_mesh (const Int np, const Int tq_order, const Int ne,
                       const Int nonunimesh, const MeshType::Enum mesh_type,
                       Mesh::GridRotation* grid_rotation, Mesh& m) {
  m.np = MeshType::is_subcell(mesh_type) ? 2 : np;
  m.tq_order = tq_order;
  m.nonuni = nonunimesh;
  if (MeshType::is_subcell(mesh_type)) {
    const auto basis = mesh_type == MeshType::gllsubcell ?
      m.basis : Basis::create(Basis::Type::uniform_offset_nodal);
    mesh::make_cubedsphere_subcell_mesh(ne, np, *basis, m.geo_p, m.geo_c2n);
  } else {
    mesh::make_cubedsphere_mesh(m.geo_p, m.geo_c2n, ne);
  }
  if (m.nonuni) mesh::make_nonuniform(m.geo_p);
  if (grid_rotation && grid_rotation->angle != 0) {
    mesh::rotate_grid(grid_rotation->axis, grid_rotation->angle,
                      grid_rotation->R, m.geo_p);
    m.grid_rotation = *grid_rotation;
  }
  mesh::make_cgll_from_geo(m.geo_p, m.geo_c2n, m.np, *m.basis, m.cgll_p, m.cgll_c2n);
  mesh::make_dgll_from_cgll(m.cgll_p, m.cgll_c2n, m.dglln2cglln, m.dgll_c2n);
  mesh::make_io_cgll_from_internal_cgll(m.cgll_p, m.cgll_c2n, m.cgll_io_c2n);
  fill_normals(m.geo_p, m.geo_c2n, m.geo_nml, m.geo_c2nml);
}

static std::string
tracer_name (const gallery::InitialCondition::Shape& shape, const Int i) {
  std::stringstream ss;
  ss << gallery::InitialCondition::to_string(shape) << i;
  return ss.str();
}

/* Two tracers, (s)rc tracer and (m)anufactured tracer, are paired. m = s at
   time 0. Until time 6, s is gradually removed from m until the exact answer is
   0. After time 6, it is gradually added until the exact answer is s.
     m(0) = s(0)
     m(t) = (1 + cos(2 pi t/T))/2 s(t)
     Ds(t)/Dt(t) = 0
     Dm(t)/Dt(t) = -pi/T sin(2 pi t/T) s(t) + (1 + cos(2 pi t/T))/2 Ds(t)/Dt
                 = -pi/T sin(2 pi t/T) s(t)
     dm(t) = s(t) int_t^{t+dt} -pi/T sin(2 pi t/T) dt
           = (cos(2 pi (t+dt)/T) - cos(2 pi t/T))/2 s(t).
   Alternatively, if we use
     m(t) = (1 - cos(2 pi t/T))/2 s(t)
     dm(t) = -(cos(2 pi (t+dt)/T) - cos(2 pi t/T))/2 s(t),
   we can start with m(0) = 0 and reach m(T/2) = s(T/2) to do a midpoint error
   calculation.
*/
class MimicSrcTerm {
  Int s_idx, m_idx, sign;
  Real fac;

public:
  typedef std::shared_ptr<MimicSrcTerm> Ptr;

  MimicSrcTerm (const Int s_idx_, const Int m_idx_, const Int sign_=1)
    : s_idx(s_idx_), m_idx(m_idx_), sign(sign_)
  {}

  Int get_s_idx () const { return s_idx; }
  Int get_m_idx () const { return m_idx; }

  void set_time (const Real& t /* sec */, const Real& dt) {
    static constexpr Real twelve_days = 3600*24*12;
    fac = sign*(std::cos(2*M_PI*(t+dt)/twelve_days) - std::cos(2*M_PI*t/twelve_days))/2;
  }

  Real get_increment (const Real& q_s) { return fac*q_s; }
};

class SrcTermMgr {
  Int cl_idx_, cl2_idx_;
  Array<Real> lat, lon;
  Array<Int> dgll2cgll;
  std::vector<MimicSrcTerm::Ptr> mimic_pairs;

  void get_toychem_idxs (const std::vector<gallery::InitialCondition::Shape>& ics) {
    cl_idx_ = cl2_idx_ = -1;
    for (size_t iic = 0; iic < ics.size(); ++iic) {
      if (tracer_name(ics[iic], 0) == "toychem10") cl_idx_ = iic;
      else if (tracer_name(ics[iic], 0) == "toychem20") cl2_idx_ = iic;
    }
    SIQK_THROW_IF(cl_idx_ != -1 && cl2_idx_ != cl_idx_+1,
                  "For convenience, require cl2_idx = cl_idx+1.");
  }

  void get_mimic_pairs (const std::vector<gallery::InitialCondition::Shape>& ics,
                        const Int cl_idx, const Int cl2_idx, const bool check_midpoint) {
    const Int nt = ics.size();
    assert(cl2_idx == cl_idx + 1);
    assert(nt - 1 - cl2_idx <= cl_idx);
    const Int sign = check_midpoint ? -1 : 1;
    for (Int i = 0; i < nt - cl2_idx - 1; ++i)
      mimic_pairs.push_back(std::make_shared<MimicSrcTerm>(i, cl2_idx + i + 1, sign));
  }

  static void get_increment (Real lat, Real lon, Real rho, Real cl, Real cl2, Real dt,
                             Real& cl_f, Real& cl2_f, const bool is_Q) {
    const auto cl_q = is_Q ? cl/rho : cl;
    const auto cl2_q = is_Q ? cl2/rho : cl2;
    slmm::gallery::ToyChem::tendency(lat, lon,
                                     cl_q, cl2_q, dt,
                                     cl_f, cl2_f);
    const auto rho_fac = is_Q ? rho : 1;
    cl_f *= rho_fac*dt;
    cl2_f *= rho_fac*dt;
  }

public:
  typedef std::shared_ptr<SrcTermMgr> Ptr;

  SrcTermMgr (const std::vector<gallery::InitialCondition::Shape>& ics,
              const AVec3s& cgll_p,
              const AIdxArray& dgll2cgll_,
              const bool check_midpoint) {
    get_toychem_idxs(ics);
    if (cl_idx_ == -1) return;
    const Int ncgll = nslices(cgll_p);
    lat.optclear_and_resize(ncgll); lon.optclear_and_resize(ncgll);
    for (Int i = 0; i < ncgll; ++i) {
      const auto n = slice(cgll_p, i);
      xyz2ll(n[0], n[1], n[2], lat[i], lon[i]);
    }
    dgll2cgll.optclear_and_resize(len(dgll2cgll_));
    for (size_t i = 0; i < dgll2cgll.size(); ++i) dgll2cgll[i] = dgll2cgll_[i];
    get_mimic_pairs(ics, cl_idx_, cl2_idx_, check_midpoint);
  }

  SrcTermMgr (const std::vector<gallery::InitialCondition::Shape>& ics,
              const Array<Real>& lat_, const Array<Real>& lon_,
              const bool check_midpoint)
    : lat(lat_), lon(lon_)
  {
    get_toychem_idxs(ics);
    if (cl_idx_ == -1) return;
    get_mimic_pairs(ics, cl_idx_, cl2_idx_, check_midpoint);
  }

  Int nsrc () const { return active() ? 2 + mimic_pairs.size() : 0; }
  Int ndof () const { return dgll2cgll.empty() ? lat.size() : dgll2cgll.size(); }
  bool active () const { return cl_idx_ >= 0 && cl2_idx_ >= 0; }
  Int tracer_idx () const { return cl_idx_; }

  void assess (const Int np2, const Real* F, const Real* rho,
               const Real* tracers /* mixing ratio */, const Int len) const {
    const auto c = slmm::gallery::ToyChem::constant;
    const Real* q1 = tracers + cl_idx_*len;
    const Real* q2 = tracers + cl2_idx_*len;
    Real li_num = 0, li_den = c, l2_num = 0, l2_den = 0, mass0 = 0, massf = 0;
    for (Int i = 0; i < len; ++i) {
      const Real val = (q1[i] + q2[i]);
      li_num = std::max(li_num, std::abs(val - c));
      l2_num += F[i]*square(val - c);
      l2_den += F[i]*square(c);
      mass0 += F[i]*c;
      massf += F[i]*rho[i]*val;
    }
    printf("toy %11.6e | %11.6e %11.6e\n",
           reldif(mass0, massf),
           std::sqrt(l2_num/l2_den), li_num/li_den);
  }

  void add_tendency (const Real& t, const Real& dt, const Real* rho,
                     Real* tracers, // mixing ratio (is_Q false) or density (is_Q true)
                     const bool is_Q) {
    const Int len = ndof();
    Real* cl = tracers + cl_idx_*len;
    Real* cl2 = tracers + cl2_idx_*len;
    assert(mimic_pairs.empty() || ! is_Q);
    for (auto& e : mimic_pairs) e->set_time(t, dt);
#   pragma omp parallel for
    for (Int i = 0; i < len; ++i) {
      const Int k = dgll2cgll.empty() ? i : dgll2cgll[i];
      Real cl_f, cl2_f;
      get_increment(lat[k], lon[k], rho[i], cl[i], cl2[i], dt, cl_f, cl2_f, is_Q);
      cl[i] += cl_f;
      cl2[i] += cl2_f;
      for (auto& e : mimic_pairs)
        tracers[e->get_m_idx()*len + i] +=
          e->get_increment(tracers[e->get_s_idx()*len + i]);
    }
  }

  void get_increment (const Real& t, const Real& dt,
                      const Real* rho, const Real* tracers, Real* tends,
                      const bool is_Q) {
    const Int len = ndof();
    const Real* cl = tracers + cl_idx_*len;
    const Real* cl2 = tracers + cl2_idx_*len;
    Real* const cl_fs = tends;
    Real* const cl2_fs = tends + len;
    assert(mimic_pairs.empty() || ! is_Q);
    for (auto& e : mimic_pairs) e->set_time(t, dt);
#   pragma omp parallel for
    for (Int i = 0; i < len; ++i) {
      const Int k = dgll2cgll.empty() ? i : dgll2cgll[i];
      get_increment(lat[k], lon[k], rho[i], cl[i], cl2[i], dt, cl_fs[i], cl2_fs[i], is_Q);
      for (auto& e : mimic_pairs)
        tends[(e->get_m_idx() - cl_idx_)*len + i] =
          e->get_increment(tracers[e->get_s_idx()*len + i]);
    }
  }
};

template <typename T> T nzvalor1 (const T& val) { return val == 0 ? 1 : val; }

class Observer {
  struct Measurements {
    Real mass, mass_sphere, min, max, mass_redistributed, mass_discrepancy;

    Measurements ()
      : mass(0), mass_sphere(0), min(1e20), max(-1),
        mass_redistributed(0), mass_discrepancy(0)
    {}
  };

  struct FinalMeasurements {
    Real l1_err, l2_err, linf_err, mass_redistributed, mass_discrepancy;

    FinalMeasurements ()
      : l1_err(0), l2_err(0), linf_err(0),
        mass_redistributed(0), mass_discrepancy(0)
    {}
  };

  struct Field {
    std::string name;
    Array<Real> ic;
    Measurements im;
    std::vector<Measurements> ms;
    FinalMeasurements fm;
  };

  const Input in_;
  const Int len_;
  std::vector<std::shared_ptr<Field> > field_;
  std::vector<Real> time_;

public:
  typedef std::shared_ptr<Observer> Ptr;
  enum Language { matlab, python };
  
  Observer (const Input& in, const Int len) : in_(in), len_(len) {}

  // Add a field. The return value is the handle of that field.
  Int add_field (const std::string& name, const Real* ic) {
    field_.push_back(std::make_shared<Field>());
    auto& f = *field_.back();
    f.name = name;
    f.ic.optclear_and_resize(len_);
    copy(f.ic.data(), ic, len_);
    return field_.size() - 1;
  }

  void set_time (const Real time) { time_.push_back(time); }

  void add_obs (const Int& field_idx, const Mesh& m,
                const ARealArray& F_gll,
                const ARealArray& F_sphere, 
                // Set q to null if you're measuring just density.
                const Real* const q, const Real* const rho,
                const Real total_mass_redistribution,
                const Real total_mass_discrepancy,
                // If final, then measure error w.r.t. ICs.
                const bool final_step = false) {
    if (field_idx >= static_cast<Int>(field_.size()))
      throw std::runtime_error("Invalid field_idx.");
    auto& f = *field_[field_idx];

    if (time_.empty())
      throw std::runtime_error("Need to set_time.");
    const bool first = time_.size() == 1;
    if (first && final_step)
      throw std::runtime_error("Both first and final_step.");
    if ( ! first) f.ms.push_back(Measurements());
    auto& meas = first ? f.im : f.ms.back();

    // Longitude mask for CBandSC.
    Array<char> lonmask;
    const bool lon_mask = (final_step &&
                           f.name.substr(0, f.name.size()-1) == "cbandsc");
    if (lon_mask) {
      const Int n = nslices(m.cgll_p);
      lonmask.optclear_and_resize(n);
#     pragma omp parallel for
      for (Int i = 0; i < n; ++i) {
        const auto n = slice(m.cgll_p, i);
        Real lat, lon;
        xyz2ll(n[0], n[1], n[2], lat, lon);
        lonmask[i] = lon >= (M_PI/6) && lon <= 5*(M_PI/6);
      }
    }

    const auto np2 = square(m.np);
    const auto err_red = [&] (const Int ci, Real* a) {
      Real meas_mass = 0, den_l1_err = 0, den_l2_err = 0,
        fm_l1_err = 0, fm_l2_err = 0, meas_mass_sphere = 0;
      for (Int j = 0, basis_idx = 0; j < m.np; ++j)
        for (Int i = 0; i < m.np; ++i, ++basis_idx) {
          const Int k = ci*np2 + basis_idx;
          const Int idx = k;
          const Real Q = (q ? q[idx] : 1)*rho[idx];
          {
            const Real w = F_gll[k];
            meas_mass += w * Q;
            if (final_step && ( ! lon_mask || lonmask[m.dglln2cglln[k]])) {
              const Real v = q ? q[idx] : rho[idx];
              const Real v0 = f.ic[idx];
              den_l1_err += w * std::abs(v0);
              den_l2_err += w * square(v0);
              fm_l1_err += w * std::abs(v - v0);
              fm_l2_err += w * square(v - v0);
            }
          }
          {
            const Real w = F_sphere[k];
            meas_mass_sphere += w * Q;
          }
        }
      a[0] += meas_mass; a[1] += den_l1_err; a[2] += den_l2_err; a[3] += fm_l1_err;
      a[4] += fm_l2_err; a[5] += meas_mass_sphere;
    };
    Real a[6];
    accumulate_threaded_bfb<6>(err_red, nslices(m.geo_c2n), a);
    const Real meas_mass = a[0], den_l1_err = a[1], den_l2_err = a[2], fm_l1_err = a[3],
      fm_l2_err = a[4], meas_mass_sphere = a[5];
    Real fm_linf_err = 0, den_linf_err = 0;
    if (final_step) {
#     pragma omp parallel for reduction (max: den_linf_err, fm_linf_err)
      for (Int ci = 0; ci < nslices(m.geo_c2n); ++ci)
        for (Int j = 0, basis_idx = 0; j < m.np; ++j)
          for (Int i = 0; i < m.np; ++i, ++basis_idx) {
            const Int k = ci*np2 + basis_idx;
            const Int idx = k;
            if (lon_mask && ! lonmask[m.dglln2cglln[k]]) continue;
            const Real v = q ? q[idx] : rho[idx];
            const Real v0 = f.ic[idx];
            den_linf_err = std::max(den_linf_err, std::abs(v0));
            fm_linf_err = std::max(fm_linf_err, std::abs(v - v0));
          }
    }

    meas.mass = meas_mass;
    meas.mass_sphere = meas_mass_sphere;
    if (final_step) {
      f.fm.l1_err = fm_l1_err / nzvalor1(den_l1_err);
      f.fm.l2_err = std::sqrt(fm_l2_err / nzvalor1(den_l2_err));
      f.fm.linf_err = fm_linf_err / nzvalor1(den_linf_err);
    }

    Real meas_min = 1, meas_max = 0;
#   pragma omp parallel for reduction (min: meas_min)
    for (Int idx = 0; idx < len_; ++idx) {
      const Real Q = q ? q[idx] : rho[idx];
      meas_min = std::min(meas_min, Q);
    }
#   pragma omp parallel for reduction (max: meas_max)
    for (Int idx = 0; idx < len_; ++idx) {
      const Real Q = q ? q[idx] : rho[idx];
      meas_max = std::max(meas_max, Q);
    }
    meas.min = meas_min;
    meas.max = meas_max;

    meas.mass_redistributed = total_mass_redistribution;
    f.fm.mass_redistributed += total_mass_redistribution;
    meas.mass_discrepancy = total_mass_discrepancy;
    f.fm.mass_discrepancy += total_mass_discrepancy;
  }

  void write (const Output& out, const std::string& fn_prefix,
              const Language& language = matlab, const bool summary = false) const
  {
    const auto fn = fn_prefix + (summary ? "_sum." : "_rit.") +
      (language == matlab ? "m" : "py");
    std::ofstream os(fn.c_str(), std::ofstream::out);
    os.precision(16);
    if ( ! os.is_open())
      throw std::runtime_error("Could not open file " + fn);

    const std::string
      row_start = language == matlab ? "" : "[",
      row_end = language == matlab ? "" : "],",
      comment = language == matlab ? "%" : "#",
      list_start = language == matlab ? "{" : "[",
      list_end = language == matlab ? "}" : "]";

    // Write header.
    if (language == matlab)
      os << "function s = ObserverOutput ()\n";
    else
      os << "def ObserverOutput():\n  import numpy as np\n  class Struct():\n"
         << "    pass\n  s = Struct()\n";
    os << "  s.ne = " << in_.ne << ";\n";
    os << "  s.np = " << (MeshType::is_subcell(in_.mesh_type) ? -in_.np : in_.np)
       << ";\n";
    os << "  s.tq_order = " << in_.tq_order << ";\n";
    os << "  s.timeint = " << TimeInt::convert(in_.timeint);
    if (TimeInt::is_interp(in_.timeint)) os << " " << in_.timeint_v_np;
    os << ";\n";
    os << "  s.nsteps = " << in_.nsteps << ";\n";
    os << "  s.ode = '" << in_.ode << "';\n";
    os << "  s.T = " << in_.T << ";\n";
    os << "  s.method = " << in_.integrate_options.method << ";\n";
    os << "  s.mono = " << in_.integrate_options.filter << ";\n";
    os << "  s.limiter = " << in_.integrate_options.limiter << ";\n";
    os << "  s.xyz_form = " << in_.xyz_form << ";\n";
    os << "  s.d2c = " << in_.integrate_options.d2c << ";\n";
    os << "  s.dmc = '" << Dmc::convert(in_.integrate_options.dmc) << "';\n";
    os << "  s.mesh = '" << MeshType::convert(in_.mesh_type) << "';\n";
    
    // Write data field names.
    os << "  s.fields = " << list_start;
    for (auto f : field_) os << " '" << f->name << "',";
    os << list_end << ";\n";
    os << "  s.measurements = " << list_start
       << "'{Q,rho}_mass', '{Q,rho}_mass_sphere', '{q,rho}_min', '{q,rho}_max'"
       << list_end << ";\n";

    // Write measurements.
    os << "  " << comment
       << " First row is data; subsequent are deviations: data(t) - data(t=0).\n";
    if (language == matlab)
      os << "  s.data = [";
    else
      os << "  s.data = np.array([";
    os << row_start;
    os << std::setw(13) << time_[0];
    for (auto f : field_)
      os << ", " << std::setw(23) << f->im.mass
         << ", " << std::setw(23) << f->im.mass_sphere
         << ", " << std::setw(23) << f->im.min
         << ", " << std::setw(23) << f->im.max;
    os << row_end << "\n";
    for (size_t ti = summary ? time_.size() - 2 : 0;
         ti < time_.size() - 1; ++ti) {
      os << row_start;
      os << std::setw(23) << time_[ti+1];
      for (auto f : field_)
        os << ", " << std::setw(23) << ((f->ms[ti].mass - f->im.mass) /
                                       f->im.mass)
           << ", " << std::setw(23) << ((f->ms[ti].mass_sphere - f->im.mass_sphere) /
                                       f->im.mass_sphere)
           << ", " << std::setw(23) << (f->ms[ti].min - f->im.min)
           << ", " << std::setw(23) << (f->ms[ti].max - f->im.max);
      os << row_end << "\n";
    }
    if (language == matlab)
      os << "];\n";
    else
      os << "]);\n";
    
    // Write final measurements.
    os << "  s.l1_err = [";
    for (auto f : field_) os << f->fm.l1_err << ", ";
    os << "];\n";
    os << "  s.l2_err = [";
    for (auto f : field_) os << f->fm.l2_err << ", ";
    os << "];\n";
    os << "  s.linf_err = [";
    for (auto f : field_) os << f->fm.linf_err << ", ";
    os << "];\n";
    os << "  s.mass_redistributed = [";
    for (auto f : field_) os << f->fm.mass_redistributed << ", ";
    os << "];\n";
    os << "  s.mass_discrepancy = [";
    for (auto f : field_) os << f->fm.mass_discrepancy << ", ";
    os << "];\n";

    if (language == python)
      os << "  s.et = Struct()\n  s.mem = Struct()\n";
    os << "  s.et.timestep = " << out.et_timestep << ";\n";
    os << "  s.mem.highwater = " << out.mem_hightwater << ";\n";
    if (language == matlab)
      os << "end\n";
    else
      os << "  return s\n";
  }

  void check (const SrcTermMgr& srcterm, const bool lauritzen_diag) const {
    bool ok = true;
    for (size_t fi = 0; fi < field_.size(); ++fi) {
      auto& f = field_[fi];
      Real mass_err = 0, min_err = 0, max_err = 0;
      const bool mimic_src_term = srcterm.active() && (Int) fi >= srcterm.tracer_idx() + 3;
      for (size_t ti = mimic_src_term ? (Int) time_.size() - 2 : 0;
           ti < time_.size() - 1; ++ti) {
        if (Dmc::use_homme_mass(in_.integrate_options.dmc))
          mass_err = std::max(
            mass_err, std::abs((f->ms[ti].mass - f->im.mass)/
                               nzvalor1(f->im.mass)));
        else
          mass_err = std::max(
            mass_err, std::abs((f->ms[ti].mass_sphere - f->im.mass_sphere)/
                               nzvalor1(f->im.mass_sphere)));
        if (f->name == "rho" ||
            Filter::is_positive_only(in_.integrate_options.filter)) {
          const Real d = f->im.min;
          if (d < 0) min_err = std::max(min_err, -d);
        } else {
          Real d = f->im.min - f->ms[ti].min;
          if (d > 0) min_err = std::max(min_err, d);
          d = f->ms[ti].max - f->im.max;
          if (d > 0) max_err = std::max(max_err, d);
        }
      }
      printf("C %3s | %1.1e | %8.1e %1.1e | %1.3e %1.3e %1.3e | %8.2e %8.2e\n",
             f->name.substr(0, 3).c_str(),
             mass_err,
             min_err, max_err,
             f->fm.l1_err, f->fm.l2_err, f->fm.linf_err,
             f->fm.mass_redistributed, f->fm.mass_discrepancy);
      if (f->name.find("toychem") == std::string::npos && ! mimic_src_term)
        ok = ok && ((in_.integrate_options.dmc == Dmc::none || mass_err < 1e-12) &&
                    ( ! in_.integrate_options.filter ||
                     std::max(min_err, max_err) < 5e-13));
    }
    std::cout << "C " << (ok ? "PASS" : "FAIL") << "\n";
    if (lauritzen_diag)
      for (auto f : field_) {
        const Real phi0 = f->ms.front().max - f->ms.front().min;
        printf("L %3s l1 %8.2e l2 %8.2e linf %8.2e phimin %9.2e phimax %9.2e\n",
               f->name.substr(0, 3).c_str(),
               f->fm.l1_err, f->fm.l2_err, f->fm.linf_err,
               (f->ms.back().min - f->ms.front().min) / phi0,
               (f->ms.back().max - f->ms.front().max) / phi0);
      }
  }

  // Hack! Don't use this for anything other than the vortex problem. That
  // problem is run w/o returning to its IC, so we need to hack our
  // infrastructure to make the "IC" the solution at the desired time.
  Real* get_ic_field_for_vortex_problem_to_overwrite (const int idx) {
    return field_[idx]->ic.data();
  }
};

static void write_binary (const Mesh& m, const AVec3s& departure_p,
                          const RemapData& rd) {
  ARealArray Js;
  calc_node_jacobians(m, departure_p, Js);
  write_binary(rd.fmm().get_M(), rd.T(), rd.dgbfi().data(),
               rd.dgbfi_gll().data(), rd.Jt().data(), Js.data(),
               m.dglln2cglln.data(), "binary.dat");
}

static void vortex_problem_hook (
  const Mesh& m, const Observer::Ptr& so, const D2Cer::Ptr& d2cer,
  const IntegrateOptions& opts, const vis::VisWriter::Ptr& iw,
  const int field_idx, const Real time)
{
  const auto dnn = len(m.dglln2cglln), cnn = nslices(m.cgll_p);
  std::vector<Real> wrk(cnn);
  for (Int i = 0; i < cnn; ++i) {
    const auto n = slice(m.cgll_p, i);
    Real lat, lon;
    xyz2ll(n[0], n[1], n[2], lat, lon);
    gallery::MovingVortices::calc_tracer(time, 1, &lat, &lon, &wrk[i]);
  }
  Real* const soln = so->get_ic_field_for_vortex_problem_to_overwrite(field_idx);
  d2cer->c2d(wrk.data(), soln);
  if (iw) iw->write(wrk.data());
}

static pg::PhysgridData::Ptr
init_physgrid (const PRefineData& pr, const Int nphys,
               const Real* rho, const Real* q,
               const int nq, const Limiter::Enum limiter) {
  const Int ncell = nslices(pr.m_t->geo_c2n), nphys2 = square(nphys);

  const auto pg = std::make_shared<pg::PhysgridData>();
  pg->ops = std::make_shared<pg::PhysgridOps>(
    *pr.m_f, pr.m_f->np, nphys, pg::Fv2Gll::Type::elrecon, pr.rd_f->Jt_.data());
  pg->limiter = limiter;
  const Int ndof = ncell*nphys2;
  pg->lat.optclear_and_resize(ndof); pg->lon.optclear_and_resize(ndof);
  pg->rho.optclear_and_resize(ndof); pg->q.optclear_and_resize(ndof*nq);

  { // Get subcell centers.
    const auto m = *pr.m_t;
#   pragma omp parallel for
    for (Int ie = 0; ie < ncell; ++ie) {
      const auto cell = slice(m.geo_c2n, ie);
      for (Int i = 0; i < nphys; ++i) {
        for (Int j = 0; j < nphys; ++j) {
          Real c[3] = {0};
          for (Int ci = 0; ci < 2; ++ci)
            for (Int cj = 0; cj < 2; ++cj) {
              Real s[3];
              siqk::sqr::calc_ref_to_sphere(
                m.geo_p, cell, 2*Real(j+cj)/nphys - 1, 2*Real(i+ci)/nphys - 1, s);
              for (Int d = 0; d < 3; ++d) c[d] += s[d];
            }
          for (Int d = 0; d < 3; ++d) c[d] *= 0.25;
          const Int os = ie*nphys2 + i*nphys + j;
          xyz2ll(c[0], c[1], c[2], pg->lat[os], pg->lon[os]);
        }
      }
    }
  }

  return pg;
}

static void integrate (
  const PRefineData::Ptr& pref, const MeshIntegrator::Ptr& mi, const Real days,
  const Int nsteps_per_12days, const std::vector<gallery::InitialCondition::Shape>& ics,
  const std::string& out_fn, const Int write_every, const IntegrateOptions& opts,
  Output& out, const Input::Debug& debug, const Observer::Ptr& so,
  const Snapshot::Ptr& snapshot)
{
  const Mesh& m = *pref->m_t;
  auto rd = pref->rd_t;
  auto md = pref->md_t;

  Timer::start(Timer::ts_setup);
  const Int dnn = len(m.dglln2cglln), cnn = nslices(m.cgll_p),
    len = dnn, nics = static_cast<Int>(ics.size());
  assert(nics > 0);

  // Initialize I/O.
  io::NetcdfWriter::Ptr ncw;
  vis::VisWriter::Ptr iw;
  if (write_every > 0) {
    if (opts.io_type == IOType::netcdf) {
      ncw = std::make_shared<io::NetcdfWriter>(
        m.cgll_p, m.cgll_io_c2n, out_fn + ".g", m.np, opts.filter);
      {
        ARealArray dgbfi, cgbfi;
        calc_gll_basis_function_integrals(m, GLL(), dgbfi, cgbfi);
        ncw->add_and_write_nodal_static_field("cgbfi_gll", cgbfi.data());
        calc_basis_function_integrals(m, m.geo_p, dgbfi, cgbfi);
        ncw->add_and_write_nodal_static_field("cgbfi_sphere", cgbfi.data());
      }
      ncw->add_nodal_field("density");
      for (Int iic = 0; iic < nics; ++iic)
        ncw->add_nodal_field(tracer_name(ics[iic], iic));
      ncw->end_definition();
    } else if (opts.io_type == IOType::internal) {
      const auto res = opts.output_resolution;
      vis::MapToArray::Ptr map;
      if (opts.vortex_problem) {
        Real xhat[] = {-1,0,0}, yhat[] = {0,0,1};
        siqk::SphereGeometry::normalize(xhat);
        map = std::make_shared<vis::BilinGLLToOrthographic>(
          m.cgll_p, m.cgll_io_c2n, xhat, yhat, 2*res+1, 2*res+1, opts.io_recon);
      } else {
        map = std::make_shared<vis::BilinGLLToLatLon>(
          m.cgll_p, m.cgll_io_c2n, 2*res+1, 4*res+1, opts.io_recon);
      }
      iw = std::make_shared<vis::VisWriter>(map, out_fn + ".bin");
    }
  }

  if (so) so->set_time(0);

  // Initialize data and workspace.
  Array<Real> tracer[2], density[2];
  Array<Real>* tracer_p[2], * density_p[2];
  for (Int i = 0; i < 2; ++i) {
    tracer[i].optclear_and_resize(nics*len);
    tracer_p[i] = &tracer[i];
    density[i].optclear_and_resize(len);
    density_p[i] = &density[i];
  }
  for (Int k = 0; k < 2; ++k)
    for (Int i = 0; i < len; ++i)
      (*density_p[k])[i] = 1;
  if (so) {
    so->add_field("rho", density_p[0]->data());
    so->add_obs(0, m, rd->dgbfi_gll(), rd->dgbfi(), nullptr,
                density_p[0]->data(), 0, 0);
  }
  Array<Real> wrk(dnn);
  // Record the initial and final states.
  Array<Real> error_data(4*len);

  const auto d2cer = std::make_shared<D2Cer>(m.dglln2cglln, rd->dgbfi_mass());

  {
    // Get the initial conditions.
    Array<Real> lat(cnn), lon(cnn);
    for (Int i = 0; i < cnn; ++i) {
      const auto n = slice(m.cgll_p, i);
      xyz2ll(n[0], n[1], n[2], lat[i], lon[i]);
    }
    if (ncw) ncw->advance_time_to(0);
    copy(error_data.data() + len, density_p[0]->data(), len);
    if (ncw || iw) {
      Real* data = wrk.data();
      d2cer->d2c(density_p[0]->data(), data);
      if (ncw) ncw->write_field("density", data);
      if (iw) iw->write(data);
    }
    for (Int iic = 0; iic < nics; ++iic) {
      Real* data = wrk.data();
      gallery::InitialCondition::init(ics[iic], nslices(m.cgll_p), lat.data(),
                                      lon.data(), data);
      // Record the ICs.
      d2cer->c2d(data, tracer_p[0]->data() + iic*len);
      if (iic == 0)
        copy(error_data.data(), tracer_p[0]->data(), len);
      if (ncw) ncw->write_field(tracer_name(ics[iic], iic), data);
      if (iw) {
        iw->write(data);
        if (opts.vortex_problem &&
            ics[iic] == gallery::InitialCondition::VortexTracer) {
          // Write this field again to match vortex_problem_hook writes.
          iw->write(data);
        }
      }
      if (so) {
        so->add_field(tracer_name(ics[iic], iic),
                      tracer_p[0]->data() + iic*len);
        so->add_obs(iic+1, m, rd->dgbfi_gll(), rd->dgbfi(),
                    tracer_p[0]->data() + iic*len, density_p[0]->data(), 0, 0);
      }
    }
  }

  LauritzenDiag::Ptr lauritzen_diag;
  if (so && opts.lauritzen_diag)
    lauritzen_diag = std::make_shared<LauritzenDiag>(
      nsteps_per_12days, len, ics, tracer_p[0]->data(), rd->dgbfi_mass().data(),
      opts.lauritzen_diag_io);

  if (Method::is_ci(opts.method)) {
    // Remap is done on density*tracer, but sometimes the tracer field doesn't
    // have the density rho in it.
    for (Int iic = 0; iic < nics; ++iic) {
      Real* const t = tracer_p[0]->data() + iic*len;
      for (Int i = 0; i < len; ++i)
        t[i] *= (*density_p[0])[i];
    }
  }

  const auto remapper = std::make_shared<Remapper>(pref, md, d2cer);
  if (opts.record_in_time && opts.filter != Filter::none)
    remapper->record_total_mass_redistribution(true);
  if (opts.perturb_rho)     remapper->set_rho_perturbation(opts.perturb_rho);
  if (opts.track_footprint) remapper->track_footprint(true);
  if (opts.subcell_bounds)  remapper->use_subcell_bounds();
  if (opts.fitext)          remapper->use_fit_extremum();

  pg::PhysgridData::Ptr pg;
  if (opts.pg > 0)
    pg = init_physgrid(*pref, opts.pg, density_p[0]->data(), tracer_p[0]->data(),
                       nics, opts.limiter);

  const auto srcterm = pg ?
    std::make_shared<SrcTermMgr>(ics, pg->lat, pg->lon, opts.check_midpoint) :
    std::make_shared<SrcTermMgr>(ics, m.cgll_p, m.dglln2cglln, opts.check_midpoint);
  Array<Real> srcterm_tendency;
  Real* srcterm_tendency_ptr = nullptr;
  Int srcterm_beg = nics, srcterm_end = nics;
  if (srcterm->active() && (pref->experiment > 1 || pg)) {
    srcterm_tendency.optclear_and_resize(srcterm->nsrc()*srcterm->ndof());
    srcterm_tendency_ptr = srcterm_tendency.data();
    srcterm_beg = srcterm->tracer_idx();
    srcterm_end = srcterm_beg + srcterm->nsrc();
    SIQK_THROW_IF(pref->experiment > 0 && opts.filter == Filter::none,
                  "Toy chemistry with p-refinement requires shape preservation.");
  }
  if (srcterm->active() && opts.check_midpoint) {
    printf("WARNING: Toy chemistry's midpoint check value is invalid.\n");
    // (because the midpoint check reference is computed using one large time
    // step, invalidating the chemistry)
  }

  AVec3s departure_p, rho_departure_p;
  resize(departure_p, mi->nnodes());
  if (opts.method == Method::isl)
    resize(rho_departure_p, nslices(m.geo_p));

  if (pg) {
    // Transfer data to physgrid consistently with p-refinement method. Do this
    // here to have proper data for the first toy-chemistry tendencies call.
    remapper->init_physgrid(density_p[0]->data(), tracer_p[0]->data(), nics, *pg);
  }

  // Time step.
  const Real T = day2sec(days);
  const Real ncycle = days/12.0;
  const Int nsteps = ncycle*nsteps_per_12days;
  const Real dt = T/nsteps;
  const Int last_step = nsteps - 1;

  ProgressBar progress_bar("integrate", last_step+1, 10);
  const auto step_t = siqk::tic();
  Timer::stop(Timer::ts_setup);

  Timer::start(Timer::ts);
  for (Int step = 0; step <= last_step; ++step) {
    Timer::start(Timer::ts_integrate);
    Real ts = dt*step;
    Real tf = step == last_step ? T : ts + dt;
    const Int cycle = 1 + step / nsteps_per_12days;

    const bool do_io = ((ncw || iw) &&
                        (((step+1) % write_every == 0 &&
                          cycle >= opts.io_start_cycle) ||
                         step == last_step));
    const bool do_Qq = (Method::is_ci(opts.method) &&
                        (do_io || so || step == last_step));

    if (srcterm->active()) {
      const Real* rho = pg ? pg->rho.data() : density_p[0]->data();
      Real* tracer = pg ? pg->q.data() : tracer_p[0]->data();
      if (pref->experiment > 1 || pg) {
        srcterm->get_increment(ts, dt, rho, tracer, srcterm_tendency_ptr,
                               Method::is_ci(opts.method));
        if (0 && pg) {
          static vis::VisWriter::Ptr iw_pg;
          if ( ! iw_pg && write_every > 0) {
            const Int ne = std::sqrt(nslices(m.geo_c2n)/6);
            iw_pg = pg::make_pg_vis(out_fn + "_srcterm", ne, pg->ops->nphys,
                                    opts.output_resolution);
          }
          if (do_io) {
            Real* const data = wrk.data();
            for (Int i = 0; i < 2; ++i) {
              const Int ndof = srcterm->ndof();
              const Real* const q = tracer + ndof*(srcterm->tracer_idx() + i);
              const Real* const dq = srcterm_tendency_ptr + i*ndof;
              for (Int i = 0; i < ndof; ++i) data[i] = q[i] + dq[i];
              iw_pg->write(data);
            }
          }
        }
      } else {
        srcterm->add_tendency(ts, dt, rho, tracer, Method::is_ci(opts.method));
        if (0) {
          static vis::VisWriter::Ptr iw_gll;
          if ( ! iw_gll && write_every > 0) {
            assert(iw);
            iw_gll = std::make_shared<vis::VisWriter>(iw->get_map(), out_fn + "_srcterm_gll.bin");
          }
          if (do_io) {
            Real* const data = wrk.data();
            for (Int i = 0; i < 2; ++i) {
              d2cer->d2c(tracer + srcterm->ndof()*(srcterm->tracer_idx() + i), data);
              iw_gll->write(data);
            }
          }
        }
      }
    }

    // Advect mesh vertices.
    switch (opts.stepping) {
    case IntegrateOptions::fwd:
      assert(Method::is_ci(opts.method));
      mi->integrate(ts, tf, departure_p);
      break;
    case IntegrateOptions::bwd:
      assert(Method::is_isl(opts.method));
      mi->integrate(tf, ts, departure_p);
      if (opts.method == Method::isl)
        mi->integrate(ts, tf, rho_departure_p);
      break;
    }
    Timer::stop(Timer::ts_integrate);

    // CDG, IR, or ISL.
    Timer::start(Timer::ts_remap);
    switch (opts.method) {
    case Method::ir:
    case Method::cdg:
      remapper->remap(departure_p,
                      tracer_p[0]->data(), tracer_p[1]->data(), nics,
                      density_p[0]->data(), density_p[1]->data());
      break;
    case Method::isl:
      if ( ! opts.rhot0)
        remapper->remap(rho_departure_p, nullptr, nullptr, 0,
                        density_p[0]->data(), density_p[1]->data());
      remapper->isl(departure_p,
                    tracer_p[0]->data(), tracer_p[1]->data(), nics,
                    density_p[0]->data(), density_p[1]->data(),
                    Filter::is_positive_only(opts.filter), false,
                    srcterm_tendency_ptr, srcterm_beg, srcterm_end, pg);
      break;
    case Method::pisl:
    case Method::pislu:
      remapper->isl(departure_p,
                    tracer_p[0]->data(), tracer_p[1]->data(), nics,
                    density_p[0]->data(), density_p[1]->data(),
                    Filter::is_positive_only(opts.filter), ! opts.rhot0,
                    srcterm_tendency_ptr, srcterm_beg, srcterm_end, pg);
      break;
    default: throw std::logic_error("Not a Method.");
    }

    if (opts.d2c)
      dss(*d2cer, Method::is_isl(opts.method), nics, len, wrk.data(),
          density_p[1]->data(), tracer_p[1]->data(),
          // For pisl* methods, isl routine takes care of rho DSS.
          ! Method::is_pisl(opts.method));

    Timer::stop(Timer::ts_remap);

    // I/O and error analysis. These are done with q, not Q, so need to convert
    // to q and then later back to Q.
    Timer::start(Timer::ts_rest);
    if (do_Qq)
      Qtoq(nics, len, density_p[1]->data(), tracer_p[1]->data());

    gdbg.write("T", rd->T());
    gdbg.write_p("geo_p", m.geo_p); gdbg.write_c2n("geo_c2n", m.geo_c2n);
    gdbg.write_p("departure_p", departure_p);
    if (step == debug.write_binary_at)
      write_binary(m, departure_p, *rd);
    
    // Netcdf I/O.
    if (do_io) {
      if (ncw) ncw->advance_time_to(tf);
      {
        Real* const data = wrk.data();
        d2cer->d2c(density_p[1]->data(), data);
        if (ncw) ncw->write_field("density", data);
        if (iw) iw->write(data);
        for (Int iic = 0; iic < nics; ++iic) {
          d2cer->d2c(tracer_p[1]->data() + iic*len, data, opts.io_no_dss);
          if (ncw) ncw->write_field(tracer_name(ics[iic], iic), data);
          if (iw) iw->write(data);
        }
      }
    }

    const bool report = (((step + 1) % nsteps_per_12days == 0) ||
                         step == last_step);

    if (opts.vortex_problem && (do_io || report)) {
      for (Int iic = 0; iic < nics; ++iic) {
        if (ics[iic] != gallery::InitialCondition::VortexTracer) continue;
        vortex_problem_hook(m, so, d2cer, opts, iw, iic+1, tf);
        break;
      }
    }

    // Observer.
    if (so) {
      if (step % nsteps_per_12days == 0)
        std::cout << "\nC cycle " << cycle << "\n";
      so->set_time(tf);
      so->add_obs(0, m, rd->dgbfi_gll(), rd->dgbfi(), nullptr,
                  density_p[1]->data(),
                  remapper->get_total_redistributed_mass(nics),
                  remapper->get_total_mass_discrepancy(nics),
                  report);
      for (Int iic = 0; iic < nics; ++iic)
        so->add_obs(iic+1, m, rd->dgbfi_gll(), rd->dgbfi(),
                    tracer_p[1]->data() + iic*len, density_p[1]->data(),
                    remapper->get_total_redistributed_mass(iic),
                    remapper->get_total_mass_discrepancy(iic),
                    report);
      if (srcterm->active())
        srcterm->assess(square(m.np), rd->dgbfi_mass().data(),
                        density_p[1]->data(), tracer_p[1]->data(), len);
      if (report) {
        printf("\n");
        so->check(*srcterm, lauritzen_diag != nullptr);
      }
      if (lauritzen_diag) {
        lauritzen_diag->run(step, tracer_p[1]->data(),
                            rd->dgbfi_mass().data(), *d2cer, out_fn);
        if (pref->experiment > 1) {
          const Real* q, * dA;
          Int len;
          remapper->get_prefine_internal(q, dA, len, step == 0);
          lauritzen_diag->run_me(len, q, dA, step == 0);
        }
        if (report) {
          lauritzen_diag->print();
          if (opts.lauritzen_diag_io)
            printf("L file %s\n", (out_fn + "-lauritzen-diag.dat").c_str());
        }
      }
    }

    if (opts.check_midpoint && snapshot &&
        // The flow direction alternates, so the midpoint snapshot is valid only
        // on every other cycle.
        step / nsteps_per_12days % 2 == 0 &&
        (step + 1) % (nsteps_per_12days/2) == 0 &&
        (step + 1) % nsteps_per_12days != 0)
      check(*snapshot, *pref->m_t, tracer_p[1]->data(), rd->dgbfi_mass().data(), dnn);

    // Record data for error analysis.
    if (step == last_step) {
      copy(error_data.data() + 2*len, tracer_p[1]->data(), len);
      copy(error_data.data() + 3*len, density_p[1]->data(), len);
    }

    if (do_Qq)
      qtoQ(nics, len, density_p[1]->data(), tracer_p[1]->data());

    swap(tracer_p);
    swap(density_p);
    progress_bar.update();
    gdbg.advance();
    gdbg.set_on(false);
  }
  const Real step_et = siqk::toc(step_t);
  Timer::stop(Timer::ts);

  const bool calc_midpoint_snapshot = ! opts.check_midpoint && snapshot != nullptr;
  out.et_timestep = step_et;
  out.mem_hightwater = siqk::get_memusage();
  if ( ! calc_midpoint_snapshot) {
    siqk::print_times("timestep", step_et);
    Timer::start(Timer::ts_error);
    const Real* const d = error_data.data();
    print_error(m, rd->dgbfi_gll(), rd->dgbfi(), d, d + len,
                d + 2*len, d + 3*len, out);
    Timer::stop(Timer::ts_error);
  }

  if (calc_midpoint_snapshot) {
    auto& s = *snapshot;
    s.m = pref->m_t;
    s.q.resize(nics);
    s.ics = ics;
    // For fields with the mimic source term, make the midpoint reference the
    // corresponding reference field.
    Int mimic_src_term_idx = -1;
    if (srcterm->active() && nics > srcterm->tracer_idx() + 2)
      mimic_src_term_idx = srcterm->tracer_idx() + 2;
    for (Int i = 0; i < nics; ++i) {
      s.q[i].resize(dnn);
      Int k = i;
      if (mimic_src_term_idx >= 0 && i >= mimic_src_term_idx) {
        k = i - mimic_src_term_idx;
        s.ics[i] = s.ics[k];
      }
      const auto d = tracer_p[0]->data() + k*dnn;
      std::copy(d, d + dnn, s.q[i].begin());
    }
  }
}

static Snapshot::Ptr make_midpoint_snapshot (
  const std::string& ode, const Int ne, const Mesh& m,
  const std::vector<gallery::InitialCondition::Shape>& ics)
{
  IntegrateOptions io;
  io.stepping = IntegrateOptions::bwd;
  io.method = Method::pislu; // 1 step with natural interpolant
  io.d2c = true;
  io.filter = Filter::none;
  io.limiter = Limiter::none;
  io.subcell_bounds = false;
  io.fitext = false;
  io.record_in_time = false;
  io.check_midpoint = false;
  io.lauritzen_diag = io.lauritzen_diag_io = false;
  io.perturb_rho = 0;
  io.dmc = Dmc::eq_homme;
  io.track_footprint = false;
  io.io_type = IOType::internal;
  io.io_no_dss = false;
  io.output_resolution = 0;
  io.pg = 0;

  const auto pref = std::make_shared<PRefineData>();
  pref->prefine = false;
  pref->experiment = 0;
  pref->m_f = pref->m_t = std::make_shared<Mesh>();
  const auto& mesh = pref->m_f;
  mesh->basis = Basis::create(Basis::Type::gll);
  //const Int plus = 4, np = m.np >= 14-plus ? 16 : m.np + plus;
  mesh->grid_rotation = m.grid_rotation;
  init_mesh(m.np, 4 /* unused */, ne, m.nonuni, MeshType::geometric,
            &mesh->grid_rotation, *mesh);
  pref->rd_f = pref->rd_t = std::make_shared<RemapData>(*mesh, io.method, io.dmc);

  const auto mi = MeshIntegratorFactory::create(ode, false, mesh->cgll_p);
  mi->config_for_best_possible_accuracy();

  Output out;
  Input::Debug dbg;
  const auto snapshot = std::make_shared<Snapshot>();
  integrate(pref, mi, 6, 2, ics, "", 0, io, out, dbg, nullptr, snapshot);

  return snapshot;
}

static Basis::Ptr parse_basis (const Method::Enum method, const std::string& basis) {
  if (Method::is_ci(method) || ! Method::is_stabilized(method))
    return Basis::create(Basis::Type::gll);
  else if (basis == "" || basis == "GllNodal")
    return Basis::create(Basis::Type::gll_nodal);
  else if (basis == "GllOffsetNodal")
    return Basis::create(Basis::Type::gll_offset_nodal);
  else if (basis == "UniformOffsetNodal")
    return Basis::create(Basis::Type::uniform_offset_nodal);
  else if (basis == "Gll")
    return Basis::create(Basis::Type::gll);
  else {
    const auto b = Basis::create_basis_from_string(basis);
    SIQK_THROW_IF( ! b, "Did not get a basis from: " << b);
    return b;
  }
}

static void run (const Input& in) {
  std::vector<gallery::InitialCondition::Shape> ics;
  for (const auto& e : in.initial_conditions)
    for (Int i = 0; i < e.n; ++i)
      ics.push_back(gallery::InitialCondition::from_string(e.name));

  const auto pref = std::make_shared<PRefineData>();
  pref->prefine = TimeInt::is_interp(in.timeint);
  pref->experiment = in.p_refine_experiment;

  const auto m = std::make_shared<Mesh>();
  {
    m->basis = parse_basis(in.integrate_options.method, in.basis);
    Mesh::GridRotation gr;
    if (in.rotate_grid) {
      if (in.integrate_options.vortex_problem) {
        // Follow Nair & Jablonowski (see comment above MovingVortices::eval for
        // full citation) and orient the grid so that the vortices move over 4
        // of the 8 cubed-sphere corners, thus probing the effects of the
        // corners. The grid corners are sometimes considered to be troublesome.
        gr.axis[0] = 1; gr.axis[1] = 0; gr.axis[2] = 0;
        // Fixed but random number a little under 1. See below.
        gr.angle = 0.97654321*M_PI/4;
      } else {
        // Rotate by a fixed but random amount to keep the center of the
        // background solid-body rotation off an element corner. An unmoving
        // solid-body rotation exactly on a collocation point excites a mild
        // instability, and since the center never moves, the fact that the
        // instability is mild doesn't help: after enough cycles, it will show
        // up. This is true for even np=3.
        gr.axis[0] = 0.11111; gr.axis[1] = -0.051515; gr.axis[2] = 1;
        gr.angle = 0.142314*M_PI;
      }
      siqk::SphereGeometry::normalize(gr.axis);
    }
    init_mesh(in.np, in.tq_order, in.ne, in.nonunimesh, in.mesh_type,
              in.rotate_grid ? &gr : nullptr, *m);
    pref->m_f = m;
  }

  MeshIntegrator::Ptr mi; {
    Mesh::Ptr m_coarse;
    if (pref->prefine) {
      m_coarse = std::make_shared<Mesh>();
      // We want the coarse basis to be stable to handle rho.
      m_coarse->basis = Basis::create(Basis::Type::gll_offset_nodal);
    }
    pref->m_v = m_coarse;
    const auto& p = (Method::is_ci(in.integrate_options.method) ?
                     m->geo_p :
                     (TimeInt::is_interp(in.timeint) ?
                      m_coarse->cgll_p :
                      m->cgll_p));
    switch (in.timeint) {
    case TimeInt::exact:
      mi = MeshIntegratorFactory::create(in.ode, in.xyz_form, p);
      if (in.ode_tight) mi->config_for_best_possible_accuracy();
      break;
    case TimeInt::line: {
      StudyTimeIntegratorOptions o;
      mi = StudyTimeIntegratorFactory::create(in.ode, o, p);
    } break;
    case TimeInt::interp:
    case TimeInt::interpline: {
      m_coarse->np = in.timeint_v_np;
      m_coarse->nonuni = m->nonuni;
      m_coarse->grid_rotation = m->grid_rotation;
      m_coarse->tq_order = m->tq_order;
      m_coarse->geo_p = m->geo_p;
      m_coarse->geo_nml = m->geo_nml;
      m_coarse->geo_c2n = m->geo_c2n;
      m_coarse->geo_c2nml = m->geo_c2nml;
      mesh::make_cgll_from_geo(m->geo_p, m->geo_c2n, in.timeint_v_np,
                               *m_coarse->basis, m_coarse->cgll_p,
                               m_coarse->cgll_c2n);
      mesh::make_dgll_from_cgll(m_coarse->cgll_p, m_coarse->cgll_c2n,
                                m_coarse->dglln2cglln, m_coarse->dgll_c2n);
      if (pref->experiment > 0)
        mesh::make_io_cgll_from_internal_cgll(m_coarse->cgll_p, m_coarse->cgll_c2n,
                                              m_coarse->cgll_io_c2n);
      const auto misub = in.timeint == TimeInt::interpline ?
        StudyTimeIntegratorFactory::create(in.ode, StudyTimeIntegratorOptions(),
                                           m_coarse->cgll_p) :
        MeshIntegratorFactory::create(in.ode, in.xyz_form, m_coarse->cgll_p);
      mi = pref->experiment == 0 ?
        std::make_shared<VelocityInterpolatorMeshIntegrator>(misub, m_coarse, m) :
        misub;
    } break;
    default:
      assert(0);
    }
  }

  const auto& opts = in.integrate_options;
  pref->rd_f = std::make_shared<RemapData>(*pref->m_f, opts.method, opts.dmc);
  if (pref->prefine) {
    pref->rd_v = std::make_shared<RemapData>(*pref->m_v, opts.method, opts.dmc);
    if (pref->experiment > 0) {
      // Both the nodal weights w and the nodal Jacobians get modified. The
      // first is according to the specific Basis type; the second is that we
      // prefer to use the interp of the v-grid Jacobians to the fine-grid ones,
      // as then a constant rho on the v grid maps to a const rho on the
      // fine grid.
      calc_pref_gll_quantities(*pref->m_v, *pref->rd_v, *pref->m_f, *pref->rd_f);
    }
  }

  if (opts.filter) {
    pref->md_f = std::make_shared<MonoData>(*pref->m_f, *pref->rd_f, opts.filter,
                                            opts.limiter);
    if (pref->prefine)
      pref->md_v = std::make_shared<MonoData>(*pref->m_v, *pref->rd_v, opts.filter,
                                              opts.limiter);
  }

  if (pref->experiment <= 1) {
    pref->m_t = pref->m_f;
    pref->rd_t = pref->rd_f;
    pref->md_t = pref->md_f;
  } else {
    pref->m_t = pref->m_v;
    pref->rd_t = pref->rd_v;
    pref->md_t = pref->md_v;
    assert(pref->m_t->nonuni == pref->m_v->nonuni);
  }
  assert(pref->m_t->nonuni == pref->m_f->nonuni);

  Output out;
  std::shared_ptr<Observer> so;
  if (in.integrate_options.record_in_time)
    so = std::make_shared<Observer>(in, len(pref->m_t->dglln2cglln));
  Snapshot::Ptr midpoint_snapshot;
  if (in.integrate_options.check_midpoint)
    midpoint_snapshot = make_midpoint_snapshot(in.ode, in.ne, *pref->m_t, ics);

  integrate(pref, mi, in.T, in.nsteps, ics, in.output_fn, in.write_every,
            in.integrate_options, out, in.debug, so, midpoint_snapshot);

  if (so) {
    so->write(out, in.output_fn, Observer::matlab);
    so->write(out, in.output_fn, Observer::matlab, true);
    so->write(out, in.output_fn, Observer::python);
    so->write(out, in.output_fn, Observer::python, true);
  }
  print_one_liner(in, out);
}

Input::Input (int argc, char** argv)
  : output_fn("tmp/out"), ode("divergent"), T(12), timeint(TimeInt::exact),
    ne(5), nonunimesh(0), nsteps(120), write_every(1), np(4), tq_order(18),
    mesh_type(MeshType::geometric), xyz_form(false)
{
  assert(eq("-foo", "-foo") && eq("--foo", "-foo") && eq("--foo", "-f", "--foo") &&
         eq("-f", "-f", "--foo"));
  program_name = argv[0];
  auto& iopts = integrate_options;
  iopts.stepping = IntegrateOptions::fwd;
  iopts.d2c = false;
  iopts.record_in_time = false;
  iopts.check_midpoint = false;
  iopts.lauritzen_diag = false;
  iopts.lauritzen_diag_io = false;
  iopts.perturb_rho = 0;
  iopts.dmc = Dmc::none;
  iopts.filter = Filter::none;
  iopts.limiter = Limiter::mn2;
  iopts.subcell_bounds = false;
  iopts.fitext = false;
  iopts.track_footprint = false;
  iopts.io_type = IOType::netcdf;
  iopts.io_recon = vis::Reconstruction::Enum::bilin;
  iopts.io_start_cycle = 0;
  iopts.io_no_dss = false;
  iopts.output_resolution = 64;
  iopts.pg = 0;
  iopts.rhot0 = false;
  iopts.vortex_problem = false;
  p_refine_experiment = 0;
  rotate_grid = false;
  ode_tight = false;
  bool tq_specified = false, method_specified = false;
  for (int i = 1; i < argc; ++i) {
    const std::string& token = argv[i];
    if (eq(token, "-o", "--output"))
      output_fn = argv[++i];
    else if (eq(token, "-T"))
      T = atof(argv[++i]);
    else if (eq(token, "-nsteps"))
      nsteps = atoi(argv[++i]);
    else if (eq(token, "-timeint"))
      timeint = TimeInt::convert(argv[++i]);
    else if (eq(token, "-prefine"))
      p_refine_experiment = atoi(argv[++i]);
    else if (eq(token, "-basis"))
      basis = argv[++i];
    else if (eq(token, "-ode"))
      ode = argv[++i];
    else if (eq(token, "-ic")) {
      initial_conditions.push_back(InitialCondition(argv[++i], 1));
      if (i+1 == argc) continue;
      if (argv[i+1][0] == '-') continue;
      auto& ic = initial_conditions.back();
      ic.n = atoi(argv[++i]);
    } else if (eq(token, "-mono", "--monotone"))
      iopts.filter = Filter::convert(argv[++i]);
    else if (eq(token, "-lim", "--limiter"))
      iopts.limiter = Limiter::convert(argv[++i]);
    else if (eq(token, "-subcellbounds"))
      iopts.subcell_bounds = true;
    else if (eq(token, "-fitext"))
      iopts.fitext = true;
    else if (eq(token, "-method", "--method")) {
      method_specified = true;
      iopts.method = Method::convert(argv[++i]);
    }
    else if (eq(token, "-np"))
      np = atoi(argv[++i]);
    else if (eq(token, "-pg"))
      iopts.pg = atoi(argv[++i]);
    else if (eq(token, "-mesh"))
      mesh_type = MeshType::convert(argv[++i]);
    else if (eq(token, "-tq")) {
      tq_order = atoi(argv[++i]);
      tq_specified = true;
    } else if (eq(token, "-ne"))
      ne = atoi(argv[++i]);
    else if (eq(token, "-nonunimesh"))
      nonunimesh = atoi(argv[++i]);
    else if (eq(token, "-we", "--write-every"))
      write_every = atoi(argv[++i]);
    else if (eq(token, "-io-start-cycle"))
      iopts.io_start_cycle = atoi(argv[++i]);
    else if (eq(token, "-io-type"))
      iopts.io_type = IOType::convert(argv[++i]);
    else if (eq(token, "-io-recon"))
      iopts.io_recon = vis::Reconstruction::convert(argv[++i]);
    else if (eq(token, "-io-nodss"))
      iopts.io_no_dss = true;
    else if (eq(token, "-res"))
      iopts.output_resolution = atoi(argv[++i]);
    else if (eq(token, "-xyz"))
      xyz_form = true;
    else if (eq(token, "-d2c"))
      iopts.d2c = true;
    else if (eq(token, "-dmc"))
      iopts.dmc = Dmc::convert(argv[++i]);
    else if (eq(token, "-rit", "--record_in_time"))
      iopts.record_in_time = true;
    else if (eq(token, "-midpoint-check"))
      iopts.check_midpoint = true;
    else if (eq(token, "-lauritzen"))
      iopts.lauritzen_diag = true;
    else if (eq(token, "-lauritzen-io"))
      iopts.lauritzen_diag_io = true;
    else if (eq(token, "-footprint", "--footprint"))
      iopts.track_footprint = true;
    else if (eq(token, "-write-binary-at"))
      debug.write_binary_at = atoi(argv[++i]);
    else if (eq(token, "-perturb-rho"))
      iopts.perturb_rho = atof(argv[++i]);
    else if (eq(token, "-rotate-grid"))
      rotate_grid = true;
    else if (eq(token, "-ode-tight"))
      ode_tight = true;
    else if (eq(token, "-rhot0"))
      iopts.rhot0 = true;
    else
      SIQK_THROW_IF(true, "Unrecognized option: " << token);
  }
  if ( ! tq_specified) {
    if (Dmc::is_facet(iopts.dmc))
      tq_order = (np-1)*4;
    else if ( ! Method::is_ir(iopts.method))
      tq_order = np >= 4 ? 20 : 14;
    else
      tq_order = np == 5 ? 20 : np == 4 ? 18 : np == 3 ? 14 : 8;
  }
  if (tq_order > 20) {
    std::cerr << "WARNING: Reducing tq_order from " << tq_order << " to 20.\n";
    tq_order = 20;
  }
  if ( ! method_specified) {
    // Default method to the one natural for the quadrature geometry.
    iopts.method = Dmc::is_facet(iopts.dmc) ?
      Method::cdg : Method::ir;
  }
  if (initial_conditions.empty())
    initial_conditions.push_back(InitialCondition("xyztrig", 1));
  if (Method::is_isl(iopts.method))
    iopts.stepping = IntegrateOptions::bwd;
  if (MeshType::is_subcell(mesh_type) &&
      Dmc::is_facet(iopts.dmc) &&
      Method::is_ir(iopts.method)) {
    std::cout << "WARNING: Switching to CDG; (QOF, IR) is not supported "
              << "for subcell mesh.\n";
    iopts.method = Method::cdg;
  }
  SIQK_THROW_IF(p_refine_experiment == 1 && ! TimeInt::is_interp(timeint),
                "For p_refine_experiment = 1, -timeint must be an interp type.");
  if (TimeInt::is_interp(timeint)) {
    if (np <= 3) timeint_v_np = np;
    else timeint_v_np = 4;
  } else
    p_refine_experiment = 0;
  if (p_refine_experiment <= 5 && p_refine_experiment > 0 && np == 4) {
    std::cout << "WARNING: For p_refine_experiment <= 5 and np = 4, we reset "
      " p_refine_experiment to 0.\n";
    p_refine_experiment = 0;
  }
  if ( ! iopts.record_in_time)
    iopts.lauritzen_diag = false;
#ifndef SLMMIR_LAURITZEN_DIAG
  std::cout << "Warning: Turning off Lauritzen diagnostics b/c not built.\n";
  iopts.lauritzen_diag = false;
#endif
  iopts.vortex_problem = (gallery::WindFieldType::from_string(ode) ==
                          gallery::WindFieldType::MovingVortices);
  print(std::cout);
  SIQK_THROW_IF(iopts.subcell_bounds && ! Method::is_isl(iopts.method),
                "-subcellbounds and not ISL is not supported.");
  SIQK_THROW_IF(iopts.limiter == Limiter::none &&
                ! (iopts.filter == Filter::none ||
                   iopts.filter == Filter::caas),
                "-lim none is valid only with -mono none or caas");
  SIQK_THROW_IF(iopts.pg > 0 &&
                ( ! Dmc::use_homme_mass(iopts.dmc) ||
                  ! Method::is_isl(iopts.method)),
                "-pg requires Homme mass and is impl'ed for ISL only");
  SIQK_THROW_IF(iopts.check_midpoint && nsteps % 2 != 0,
                "nsteps must be even to check the midpoint.");
  SIQK_THROW_IF(iopts.rhot0 && ! Method::is_isl(iopts.method),
                "-rhot0 is impl'ed only for ISL methods.");
}

void Input::print (std::ostream& os) const {
  os << "output (-o): " << output_fn << "\n"
     << "ode (-ode, " << MeshIntegratorFactory::get_inputs() << "): "
     <<   ode << "\n"
     << "xyz_form (-xyz): " << xyz_form << "\n"
     << "T (-T): " << T << " [day]\n"
     << "timeint " << TimeInt::convert(timeint);
  if (TimeInt::is_interp(timeint)) os << " " << timeint_v_np;
  os << "\n";
  os << "nsteps (-nsteps): " << nsteps << "\n"
     << "np (-np): " << np << "\n"
     << "pg (-pg): " << integrate_options.pg << "\n"
     << "tq (-tq): " << tq_order << "\n"
     << "ne (-ne): " << ne << "\n"
     << "nonunimesh (-nonunimesh): " << nonunimesh << "\n"
     << "method (-method): " << Method::convert(integrate_options.method) << "\n"
     << "mesh (-mesh): " << MeshType::convert(mesh_type) << "\n"
     << "dmc (-dmc, " << Dmc::get_inputs() << "): "
     <<   Dmc::convert(integrate_options.dmc) << "\n"
     << "monotone_type (-mono, {...}): "
     <<   Filter::convert(integrate_options.filter) << "\n"
     << "limiter_type (-lim, {...}): "
     <<   Limiter::convert(integrate_options.limiter) << "\n"
     << "-subcellbounds: " << integrate_options.subcell_bounds << "\n"
     << "write every (-we): " << write_every << "\n"
     << "d2c (-d2c): " << integrate_options.d2c << "\n"
     << "record_in_time (-rit): " << integrate_options.record_in_time << "\n"
     << "p-refine exp (-prefine): " << p_refine_experiment << "\n"
     << "rotate grid (-rotate-grid): " << rotate_grid << "\n"
     << "vortex problem: " << integrate_options.vortex_problem << "\n"
     << "ODE tightest tol: " << ode_tight << "\n";
  if (integrate_options.fitext)
    os << "fitext\n";
  for (auto& ic : initial_conditions)
    os << "initial condition (-ic, "
       << gallery::InitialCondition::get_inputs() << "): "
       << ic.name << " " << ic.n << "\n";
}

int main (int argc, char** argv) {
  Kokkos::initialize(argc, argv); {
    Timer::init();
    Timer::start(Timer::total);
    Input in(argc, argv);
    run(in);
    Timer::stop(Timer::total);
    Timer::print();
  } Kokkos::finalize_all();
}
