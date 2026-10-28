#ifndef INCLUDE_MPAS_CI_UTIL_HPP
#define INCLUDE_MPAS_CI_UTIL_HPP

// Cell integration utilities.

#include "mpas_mesh.hpp"
#include "mpas_siqp.hpp"
#include "mpas_comm.hpp"
#include "mpas_util.hpp"
#include "mpas_triquad.hpp"

struct TriangleQuadrature;

struct PolyClipper {
  static constexpr int max_nvert = 20;
  static constexpr int max_src_overlap = 20;

  // In and out vertex lists.
  Ar2 vi, vo;
  int ni, no;

  PolyClipper () {
    vi.reset(3, max_nvert, rbuf);
    vo.reset(3, max_nvert, rbuf + 3*max_nvert);
  }

  void clip (const siqp::sh::Mesh& m, const Int sc) {
#ifndef NDEBUG
    const bool out =
#endif
      siqp::sh::clip_against_poly(m, sc, vi, ni, vo, no, rwrk, max_nvert);
    assert(out);
  }

private:
  Real rbuf[6*max_nvert], rwrk[3*PolyClipper::max_nvert];
};

struct OTSearchFunctor : public siqp::Octree::Functor, HitsAccumulator<256> {
  bool ok;
  OTSearchFunctor () : ok(true) {}
  void operator() (const int i) override {
    const bool stat = HitsAccumulator<max_nhit>::operator()(i);
    if ( ! stat) ok = false;
  }
};

struct OctreeMgr {
  typedef std::shared_ptr<OctreeMgr> Ptr;
  typedef std::shared_ptr<const OctreeMgr> ConstPtr;
  
  siqp::Octree::ConstPtr ot;

  OctreeMgr(const Acr2& p, const Aci2& e);

  // Build bounding box for element e and run Octree::apply. Does not call
  // f.post().
  void apply(const Acr3& p_advect, const int lvl, const int* e, const int e_nmax,
             OTSearchFunctor& f) const;
};

struct CellIntegrator {
  struct Integrand {
    // (x,y,z) are unit spherical coordinates.
    virtual Real operator()(const Real x, const Real y, const Real z) const = 0;
  };

  CellIntegrator(const MpasMesh::ConstPtr& mm, const Int tq_order);
  void integrate(const Int ci, const Integrand& f, Real& area, Real& mass);

private:
  MpasMesh::ConstPtr mm;
  int tq_order;
  TriangleQuadrature tq;
};

struct CellShapeFunction {
  static constexpr bool scale_by_arc_length = false;
  static constexpr int  max_ho_ncoef = 6;
  static constexpr bool two_halo = false;
  static constexpr int  halo_edge = two_halo ? 2 : 1;
  static constexpr int  max_halo_slots = two_halo ? max_2halo_slots : max_1halo_slots;
  static constexpr bool A_integrated = false;

  bool modal;
  int order, tq_order, ho_ncoef, ncell;
  Ar2 ho_coefs, mass_consts;

  void integrate_basis_over_polygon(
    const TriangleQuadrature& tq, const Array2D<const Real>& vo, const int no,
    const Int sc, Real* basis) const;
  Real integrate_over_polygon(
    const TriangleQuadrature& tq, const Array2D<const Real>& vo, const int no,
    const Real* const shape, const Int sc) const;

  // The convention is to exclude the cell-center value, r0, from the integrals
  // giving the shape coefficients and then to add it back at the end.
  inline Real integrate (const Real r0, const Real* const shape,
                         const Real* const basis) const {
    if (modal) {
      // Coefficients are for up to x, y, x^2, x y, y^2.
      Real accum = r0*basis[0];
      for (Int i = 0; i < ho_ncoef; ++i) accum += shape[i]*basis[i+1];
      return accum;
    } else {
      // Coefficients are for up to 1, x, y, x^2, x y, y^2. For consistency in
      // other code, we have (r0 + shape[0]) associated with basis[0] here.
      Real accum = r0*basis[0];
      for (Int i = 0; i < ho_ncoef; ++i) accum += shape[i]*basis[i];
      return accum;
    }
  }

  static inline Real
  eval_shape_fn (const int order, const Real* shape, const Real* mass_const,
                 const Real x, const Real y) {
    if (mass_const) {
      if (order == 1) return 0;
      const Real f = shape[0]*(x - mass_const[0]) + shape[1]*(y - mass_const[1]);
      if (order == 2) return f;
      return (f +
              (shape[2]*(x*x - mass_const[2]) +
               shape[3]*(x*y - mass_const[3]) +
               shape[4]*(y*y - mass_const[4])));
    } else {
      Real f = shape[0];
      if (order == 1) return f;
      f += shape[1]*x + shape[2]*y;
      if (order == 2) return f;
      return (f +
              (shape[3]*(x*x) +
               shape[4]*(x*y) +
               shape[5]*(y*y)));
    }
  }

protected:
  MpasMesh::ConstPtr mm;
  Comm::Ptr he;

  void init(const MpasMesh::ConstPtr& mm, const siqp::sh::Mesh& m,
            const Comm::Ptr& he,
            const int order,      // basis function order
            const int tq_order,   // triangle quadrature order
            const bool only_own); // ncells to init, t.ncown (true) or t.ncell (false)

  void init_mass_consts(const Array2D<const Real>& p, const Array2D<const int>& e,
                        const int ci);
};

struct CellFromCellsModalShapeFunction : public CellShapeFunction {
  bool using_fv_limiter;
  Ai1 cell_nstencil;
  Ai2 c2c;

  void init(const MpasMesh::ConstPtr& mm, const siqp::sh::Mesh& m,
            const Comm::Ptr& he, const int order, const int tq_order,
            const bool using_fv_limiter, const bool only_own = false);
  void reinit(const siqp::sh::Mesh& m);

  inline void
  calc_shape (const int ci, const Acr3& f_cc, // nq,nlvl,ncell
              const int fi, const int li, Real* const shape) const {
    if (order == 1) return;
    const int* const stencil = c2c(ci);
    const int stensz = cell_nstencil[ci];
    const auto ho_coef_s = ho_coefs(ci) + 9;
    assert(stensz >= 3 && stensz <= max_halo_slots);
    for (int d = 0; d < ho_ncoef; ++d) {
      const Real* const ho_coef = ho_coef_s + stensz*d;
      Real accum = 0;
      const auto f_cc_ci = f_cc(fi,li,ci);
      for (int j = 0; j < stensz; ++j)
        accum += ho_coef[j]*(f_cc(fi,li,stencil[j]) - f_cc_ci);
      shape[d] = accum;
    }
  }

  inline void
  eval_shape_at_vtxs (const int ci, const Real& f_cc_ci, const Real* const shape,
                      const int fi, const int li, Real* const f_vtx) const {
    const auto& t = mm->t;
    const auto& vtxs = mm->v;
    const auto ne = t.c2ne[ci];
    if (order == 1) {
      for (Int i = 0; i < ne; ++i) f_vtx[i] = f_cc_ci;
      return;
    }
    const Real* const cc = ho_coefs(ci);
    const Real* const xaxis = cc + 3;
    const Real* const yaxis = cc + 6;
    for (Int i = 0; i < ne; ++i) {
      const auto vi = t.get_c2v(ci,i);
      Vec vtx; fill(vtxs, vi, vtx, true);
      Real x, y;
      calc_xy(cc, vtx, xaxis, yaxis, scale_by_arc_length, x, y);
      f_vtx[i] = f_cc_ci + eval_shape_fn(order, shape, mass_consts(ci), x, y);
    }
  }
};

struct CellFromVtxsShapeFunctionBase : public CellShapeFunction {
  inline void calc_shape (const int ci, const Real& f_cc_ci,
                          const Acr4& f_vtx, // vtx,nq,nlvl,ncell
                          const int fi, const int li, Real* const shape) const {
    const auto& t = mm->t;
    const auto ne = t.c2ne[ci];
    if (order == 1) {
      Real accum = 0;
      for (Int j = 0; j < ne; ++j)
        accum += f_vtx(j,fi,li,ci);
      *shape = accum/ne - f_cc_ci;
      return;
    }
    const auto ho_coef_s = ho_coefs(ci) + 9;
    for (Int d = 0; d < ho_ncoef; ++d) {
      const Real* const ho_coef = ho_coef_s + ne*d;
      Real accum = 0;
      for (Int j = 0; j < ne; ++j)
        accum += ho_coef[j]*(f_vtx(j,fi,li,ci) - f_cc_ci);
      shape[d] = accum;
    }
  }  
};

// Also uses the target cell's cell-center value, i.e., the cell mean.
struct CellFromCellAndVtxsModalShapeFunction : public CellFromVtxsShapeFunctionBase {
  void init(const MpasMesh::ConstPtr& mm, const siqp::sh::Mesh& m,
            const Comm::Ptr& he, const int order, const int tq_order,
            const bool only_own = false);
  void reinit(const siqp::sh::Mesh& m);
};

struct CellFromVtxsShapeFunction : public CellFromVtxsShapeFunctionBase {
  // For Eulerian grid.
  void init(const MpasMesh::ConstPtr& mm, const Comm::Ptr& he, const int order,
            const int tq_order, const bool only_own = false);
  void reinit();
  // For moving grid. Unlike the other CellShapeFunction classes, this class
  // permits cell vertices to differ from the values in mm->v. It's fine to call
  // the Eulerian-grid init at initialization and then the Lagrangian-grid
  // reinit during time stepping.
  void reinit(const Array3D<const Real>& p /* (3,nlvl,nvtx) */, const int lvl);

  // Eval the shape function at unit point p.
  Real eval_shape_fn (const int ci, const Real* shape, const Vec p) {
    const Real* const cc = ho_coefs(ci);
    const Real* const xaxis = cc + 3;
    const Real* const yaxis = cc + 6;
    Real x, y;
    calc_xy(cc, p, xaxis, yaxis, scale_by_arc_length, x, y);
    return CellShapeFunction::eval_shape_fn(order, shape, nullptr, x, y);
  }

private:
  void init0(const MpasMesh::ConstPtr& mm, const Comm::Ptr& he, const int order,
             const int tq_order, const bool only_own);
};

bool unittest_CellShapeFunction(const MpasMesh::ConstPtr& mm,
                                const siqp::sh::Mesh& m, const Comm::Ptr& he);

inline bool at_halo_edge (const Comm& he, const Int sc) {
  if ( ! he.parallel()) return false;
  // Reconstruct in cells in the halo except at the halo edge.
  const Int dist = he.cdhalo()(sc);
  return dist && dist <= CellShapeFunction::halo_edge;
}

// Some utils for testing cell-integrated methods.
namespace citest {

template <typename Fn>
struct TestIntegrand : public CellIntegrator::Integrand {
  Fn f;
  TestIntegrand (const Fn& f_) : f(f_) {}
  Real operator() (const Real x, const Real y, const Real z) const override
  { return f(x, y, z); }
};

template <typename Fn>
static TestIntegrand<Fn> create_TestIntegrand (const Fn& f) {
  return TestIntegrand<Fn>(f);
}

Real test_fn(const Int f, const Real x, const Real y, const Real z);

// Fill cell data. Comm is not needed after this call because we fill all cells.
void fill_cc_data(const MpasMesh::ConstPtr& mm, const Int fac, const Ar2& f,
                  const int tq_order = 8);

} // namespace citest

#endif
