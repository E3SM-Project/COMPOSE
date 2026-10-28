#include "mpas_ci_util.hpp"
#include "mpas_triquad.hpp"
#include "mpas_util.hpp"
#include "mpas_geometry.hpp"

OctreeMgr::OctreeMgr (const Acr2& p, const Aci2& e) {
  siqp::Octree::Options o;
  o.max_depth = 32;
  ot = std::make_shared<siqp::Octree>(p, e, o);
}

template <typename T> struct Sub02 {
  Sub02 (const Array3D<T>& a_, const int idx_omit_) : a(a_), idx_omit(idx_omit_) {}
  const T& operator() (const int r, const int c) const { return a(r,idx_omit,c); }
private:
  const Array3D<T>& a;
  Int idx_omit;
};

void OctreeMgr
::apply (const Acr3& p_advect, const int lvl, const int* e, const int e_nmax,
         OTSearchFunctor& f) const {
  assert(f.nhit == 0 && f.ok);
  // Clipped element bounding box.
  Real tbb[6];
  Sub02<const Real> p2(p_advect, lvl);
  siqp::Octree::calc_bb(p2, e, e_nmax, tbb);
  // Get list of possible overlaps.
  ot->apply(tbb, f);
}

CellIntegrator
::CellIntegrator (const MpasMesh::ConstPtr& mm_, const Int tq_order_)
  : mm(mm_), tq_order(tq_order_)
{}

void CellIntegrator
::integrate (const Int ci, const Integrand& f, Real& area, Real& mass) {
  const Real* tq_alpha, * tq_w;
  Int tqn;
  tq.get_coef(4, tq_alpha, tq_w, tqn);
  const auto& t = mm->t;
  const auto& vs = mm->v;
  const Int ne = t.c2ne[ci];
  Vec v0; fill(mm->c, ci, v0, true);
  area = 0; mass = 0;
  for (Int i = 0; i < ne; ++i) {
    const Int ie = t.get_c2e(ci,i), vi1 = t.get_e2v(ie,0), vi2 = t.get_e2v(ie,1);
    Vec v1; fill(vs, vi1, v1, true);
    Vec v2; fill(vs, vi2, v2, true);
    for (Int qi = 0; qi < tqn; ++qi) {
      Real u[3];
      const Real jac = geo::calc_tri_jacobian(v0, v1, v2, slice(tq_alpha, qi), u);
      const Real y = f(u[0], u[1], u[2]);
      area += tq_w[qi]*jac;
      mass += tq_w[qi]*jac*y;
    }
  }
  area *= 0.5; mass *= 0.5;
}

typedef CellShapeFunction sf;

static Real
integrate_over_polygon (
  const int order, const TriangleQuadrature& tq, const int tq_order,
  const Array2D<const Real>& vo, const int no,
  const Real* const mass_consts, const Real* const sc_cc,
  const Real* const xaxis, const Real* const yaxis, const Real* const shape)
{
  assert(order > 1);
  const Real* coord;
  const Real* weight;
  int nqp;
  tq.get_coef(tq_order, coord, weight, nqp);
  Real Qaccum = 0;
  for (Int i = 1; i < no-1; ++i) {
    for (Int k = 0; k < nqp; ++k) {
      Vec p;
      const Real jac = geo::calc_tri_jacobian(
        slice(vo,0), slice(vo,i), slice(vo,i+1),
        slice(coord,k), p, false);
      Real x, y;
      calc_xy(sc_cc, p, xaxis, yaxis, sf::scale_by_arc_length, x, y);
      const Real Qs_reconstruct = sf::eval_shape_fn(order, shape, mass_consts, x, y);
      Qaccum += Qs_reconstruct*weight[k]*jac;
    }
  }
  Qaccum *= 0.5;
  return Qaccum;
}

static void
lsinv_shape_function (const int order, const int n, const Real* mass_const,
                      const Real* xs, const Real* ys, Real* B) {
  Real A[sf::max_ho_ncoef*sf::max_halo_slots];
  int ncol;
  if (mass_const) { // modal
    ncol = order == 2 ? 2 : 5;
    for (int i = 0; i < n; ++i) {
      A[n*0+i] = xs[i] - mass_const[0];
      A[n*1+i] = ys[i] - mass_const[1];
      if (order == 2) continue;
      A[n*2+i] = square(xs[i]) - mass_const[2];
      A[n*3+i] = xs[i]*ys[i]   - mass_const[3];
      A[n*4+i] = square(ys[i]) - mass_const[4];
    }
  } else {
    ncol = order == 2 ? 3 : order == 3 ? 6 : 1;
    for (int i = 0; i < n; ++i) {
      A[n*0+i] = 1;
      if (order == 1) continue;
      A[n*1+i] = xs[i];
      A[n*2+i] = ys[i];
      if (order == 2) continue;
      A[n*3+i] = square(xs[i]);
      A[n*4+i] = xs[i]*ys[i];
      A[n*5+i] = square(ys[i]);
    }
  }
  Real tau[sf::max_ho_ncoef*sf::max_halo_slots], wrk[sf::max_halo_slots],
    Q[sf::max_ho_ncoef*sf::max_halo_slots];
  lsinv(A, n, ncol, tau, wrk, Q, B);
}

void CellShapeFunction
::init (const MpasMesh::ConstPtr& mm_, const siqp::sh::Mesh& m,
        const Comm::Ptr& he_, const int order_, const int tq_order_,
        const bool only_own) {
  mm = mm_;
  he = he_;
  order = order_;
  tq_order = tq_order_;
  ho_ncoef = (order == 2 ? 2 : order == 3 ? 5 : 0);
  ncell = only_own ? mm->ncown : mm->t.ncell;
  if (order < 2) return;
  mass_consts.reset(ho_ncoef, ncell);
}

// Make linear and quadratic terms each integrate to 0. The centroidal property
// makes this constant 0 for the linear terms and an exact centroid, but we find
// that it's hard to get mass conservation to machine precision without
// compensating the linear term by a little.
void CellShapeFunction
::init_mass_consts (const Array2D<const Real>& p, const Array2D<const int>& e,
                    const int ci) {
  Real* const mass_const = mass_consts(ci);
  const int ne = mm->t.c2ne[ci];
  const Real* const cc = ho_coefs(ci);
  const Real* const xaxis = cc + 3;
  const Real* const yaxis = cc + 6;
  TriangleQuadrature tq;
  Real buf[3*max_1halo_slots];
  Ar2 vo(3, ne, buf);
  for (int i = 0; i < ne; ++i) geo::copy(vo(i), p(e(i,ci)));
  Real shape[max_ho_ncoef];
  for (int j = 0; j < ho_ncoef; ++j) shape[j] = 0;
  for (int j = 0; j < ho_ncoef; ++j) mass_const[j] = 0;
  for (int j = 0; j < ho_ncoef; ++j) {
    shape[j] = 1;
    mass_const[j] = 0;
    mass_const[j] = (::integrate_over_polygon(order, tq, tq_order, vo, ne,
                                              mass_consts(ci), cc, xaxis, yaxis,
                                              shape) /
                     mm->g.my_areaCell[ci]);
    shape[j] = 0;
  }
}

void CellShapeFunction::integrate_basis_over_polygon (
  const TriangleQuadrature& tq, const Array2D<const Real>& vo, const int no,
  const Int sc, Real* basis) const
{
  if (order == 1) {
    basis[0] = geo::calc_area(tq, vo.data(), no);
    return;
  }
  const Real* const mass_const = modal ? mass_consts(sc) : nullptr;
  const Real* const B = ho_coefs(sc);
  const Real* const sc_cc = B;
  const Real* const xaxis = B + 3;
  const Real* const yaxis = B + 6;
  const Real* coord;
  const Real* weight;
  int nqp;
  tq.get_coef(tq_order, coord, weight, nqp);
  for (int i = 0; i < ho_ncoef+1; ++i) basis[i] = 0;
  for (Int i = 1; i < no-1; ++i) {
    for (Int k = 0; k < nqp; ++k) {
      Vec p;
      const Real jac = geo::calc_tri_jacobian(
        slice(vo,0), slice(vo,i), slice(vo,i+1),
        slice(coord, k), p, false);
      Real x, y;
      calc_xy(sc_cc, p, xaxis, yaxis, sf::scale_by_arc_length, x, y);
      const Real c = 0.5*weight[k]*jac;
      basis[0] += c;
      if (order == 1) continue;
      if (modal) {
        basis[1] += c*(x - mass_const[0]);
        basis[2] += c*(y - mass_const[1]);
        if (order == 2) continue;
        basis[3] += c*(x*x - mass_const[2]);
        basis[4] += c*(x*y - mass_const[3]);
        basis[5] += c*(y*y - mass_const[4]);
      } else {
        basis[1] += c*x;
        basis[2] += c*y;
        if (order == 2) continue;
        basis[3] += c*(x*x);
        basis[4] += c*(x*y);
        basis[5] += c*(y*y);
      }
    }
  }
}

Real CellShapeFunction::integrate_over_polygon (
  const TriangleQuadrature& tq, const Array2D<const Real>& vo, const int no,
  const Real* const shape, const Int sc) const
{
  if (order == 1) return shape[0]*geo::calc_area(tq, vo, no, tq_order);
  const Real* cc = ho_coefs(sc);
  return ::integrate_over_polygon(order, tq, tq_order, vo, no,
                                  modal ? mass_consts(sc) : nullptr,
                                  cc, cc+3, cc+6, shape);
}

void CellFromCellsModalShapeFunction
::init (const MpasMesh::ConstPtr& mm_, const siqp::sh::Mesh& m,
        const Comm::Ptr& he_, const int order_, const int tq_order_,
        const bool using_fv_limiter_, const bool only_own) {
  CellShapeFunction::init(mm_, m, he_, order_, tq_order_, only_own);
  modal = true;
  using_fv_limiter = using_fv_limiter_;
  if (order < 2) return;

  const int max_stencil_size = two_halo ? max_2halo_slots : max_vtx;
  c2c.reset(max_stencil_size, ncell);
  cell_nstencil.reset(ncell);
  assert(ho_ncoef > 0);
  ho_coefs.reset(3*3 + ho_ncoef*(max_halo_slots-1), ncell);

  reinit(m);
}

static void
lsinv_modal_shape_function (
  const MpasMesh& mm, const Array2D<const Real>& p, const Array2D<const int>& e,
  const int order, const int tq_order, const Real* const mass_consts,
  const Real* sc_cc, const Real* xaxis, const Real* yaxis, const int n,
  const int* stencil, const int ho_ncoef, Real* B)
{
  Real A[sf::max_ho_ncoef*sf::max_halo_slots];
  TriangleQuadrature tq;
  Real buf[3*sf::max_halo_slots];

  for (int cidx = 0; cidx < n; ++cidx) {
    const int ci = stencil[cidx];
    const auto elem = e(ci);
    int ine = 0;
    for (int i = 0; i < max_vtx; ++i) {
      if (elem[i] < 0) break;
      ++ine;
    }
    Ar2 vo(3, ine, buf);
    for (int i = 0; i < ine; ++i) geo::copy(vo(i), p(elem[i]));

    Real shape[sf::max_ho_ncoef] = {0};
    for (int j = 0; j < ho_ncoef; ++j) {
      shape[j] = 1;
      A[n*j + cidx] = (integrate_over_polygon(order, tq, tq_order, vo, ine,
                                              mass_consts, sc_cc, xaxis, yaxis,
                                              shape) /
                       mm.g.my_areaCell[ci]);
      shape[j] = 0;
    }
  }

  Real tau[sf::max_ho_ncoef*sf::max_halo_slots], wrk[sf::max_halo_slots],
    Q[sf::max_ho_ncoef*sf::max_halo_slots];
  lsinv(A, n, ho_ncoef, tau, wrk, Q, B);
}

void CellFromCellsModalShapeFunction
::reinit (const siqp::sh::Mesh& m) {
  const auto& dist = he->cdhalo();
  ompparfor for (int ci = 0; ci < ncell; ++ci) {
    if (he->parallel() && dist(ci) && (Int) dist(ci) <= halo_edge)
      continue;
    
    // Local x-y coord system.
    Vec xaxis, yaxis;
    Real xs[max_halo_slots], ys[max_halo_slots];
    int* stencil = c2c(ci);
    cell_nstencil[ci] = map_cells_to_cell_plane(*mm, ci, xaxis, yaxis, xs, ys,
                                                stencil, scale_by_arc_length,
                                                two_halo);

    Real* B = ho_coefs(ci);
    const auto cc = slice(mm->g.my_cellCentroid, ci);
    geo::copy(B, cc); geo::copy(B+3, xaxis); geo::copy(B+6, yaxis);

    Real* const mass_const = mass_consts(ci);
    if (using_fv_limiter)
      for (int j = 0; j < ho_ncoef; ++j)
        mass_const[j] = 0;
    else
      init_mass_consts(m.p, m.e, ci);

    if ( ! A_integrated)
      lsinv_shape_function(order, cell_nstencil[ci], mass_const, xs+1, ys+1, B+9);
  }

  // Least-squares operator to reconstruct linear shape function. Need complete
  // mass_const before running the integration-based one b/c (i)
  // integrate_over_polygon needs it and (ii) it's called on the cell's nbrs.
  if (A_integrated) {
    ompparfor for (int ci = 0; ci < ncell; ++ci) {
      const int* stencil = c2c(ci);
      Real* cc = ho_coefs(ci);
      const Real* xaxis = cc + 3;
      const Real* yaxis = cc + 6;
      const int n = cell_nstencil[ci];
      require((order == 2 && n >= 2) || (order == 3 && n >= 5));
      lsinv_modal_shape_function(*mm, m.p, m.e, order, tq_order,
                                 mass_consts(ci), cc, xaxis, yaxis, n, stencil,
                                 ho_ncoef, cc+9);
    }
  }
}

void CellFromCellAndVtxsModalShapeFunction
::init (const MpasMesh::ConstPtr& mm_, const siqp::sh::Mesh& m,
        const Comm::Ptr& he_, const int order_, const int tq_order_,
        const bool only_own) {
  CellShapeFunction::init(mm_, m, he_, order_, tq_order_, only_own);
  modal = true;
  if (order < 2) return;
  ho_coefs.reset(3*3 + ho_ncoef*(max_vtx-1), ncell);
  reinit(m);
}

void CellFromCellAndVtxsModalShapeFunction
::reinit (const siqp::sh::Mesh& m) {
  const auto& t = mm->t;
  ompparfor for (int ci = 0; ci < ncell; ++ci) {
    // Local x-y coord system.
    Vec xaxis, yaxis;
    Real xs[max_halo_slots], ys[max_halo_slots];
    map_vtxs_to_cell_plane(*mm, ci, xaxis, yaxis, xs, ys,
                           nullptr, scale_by_arc_length);

    Real* B = ho_coefs(ci);
    const auto cc = slice(mm->g.my_cellCentroid, ci);
    geo::copy(B, cc); geo::copy(B+3, xaxis); geo::copy(B+6, yaxis);

    Real* const mass_const = mass_consts(ci);
    init_mass_consts(m.p, m.e, ci);

    lsinv_shape_function(order, t.c2ne[ci], mass_const, xs, ys, B+9);
  }
}

void CellFromVtxsShapeFunction
::init0 (const MpasMesh::ConstPtr& mm_, const Comm::Ptr& he_, const int order_,
         const int tq_order_, const bool only_own) {
  mm = mm_;
  he = he_;
  order = order_;
  tq_order = tq_order_;
  ho_ncoef = (order == 2 ? 3 : order == 3 ? 6 : 1);
  assert(ho_ncoef <= max_ho_ncoef);
  modal = false;
  ncell = only_own ? mm->ncown : mm->t.ncell;
  if (order < 2) return;
  ho_coefs.reset(3*3 + ho_ncoef*(max_vtx-1), ncell);  
}

void CellFromVtxsShapeFunction
::init (const MpasMesh::ConstPtr& mm, const Comm::Ptr& he, const int order,
        const int tq_order, const bool only_own) {
  init0(mm, he, order, tq_order, only_own);
  if (order < 2) return;
  reinit();
}

// Get a cell-center value from vtx values by linear LS. Then find the quadratic
// shape fn's coefs using this augmented list of values.
static void
lsinv_linvtx2quadratic_shape_function (const int n, const Real* xs,
                                       const Real* ys, Real* B) {
  lsinv_shape_function(2, n, nullptr, xs, ys, B);
  // First row of this matrix gives the constant in the linear shape function,
  // which is the value at the cell center.
  Real a[max_vtx];
  for (int i = 0; i < n; ++i) a[i] = B[i];
  // Assemble vtxs and cc in the local coord system.
  Real xs1[max_vtx+1], ys1[max_vtx+1];
  for (int i = 0; i < n; ++i) xs1[i] = xs[i];
  for (int i = 0; i < n; ++i) ys1[i] = ys[i];
  xs1[n] = 0; ys1[n] = 0;
  // Quadratic shape fn from (vtx,cc).
  Real B1[6*max_vtx];
  lsinv_shape_function(3, n+1, nullptr, xs1, ys1, B1);
  // B = B1 [I; a]
  for (int row = 0; row < 6; ++row)
    for (int col = 0; col < n; ++col)
      B[n*row + col] = B1[(n+1)*row + col] + a[col]*B1[(n+1)*row + n];
}

static void
ccvtx_lsinv_shape_function (int& order, const int n,
                            const Real* xs, const Real* ys, Real* B) {
  if (order == 3) {
    // Not really order 3, but uses quadratic shape fn.
    lsinv_linvtx2quadratic_shape_function(n, xs, ys, B);
  } else
    lsinv_shape_function(order, n, nullptr, xs, ys, B);
}

void CellFromVtxsShapeFunction::reinit () {
  if (order < 2) return;
  const auto& t = mm->t;
  ompparfor for (int ci = 0; ci < ncell; ++ci) {
    // Local x-y coord system.
    Vec xaxis, yaxis;
    Real xs[max_halo_slots], ys[max_halo_slots];
    map_vtxs_to_cell_plane(*mm, ci, xaxis, yaxis, xs, ys,
                           nullptr, scale_by_arc_length);

    Real* B = ho_coefs(ci);
    const auto cc = slice(mm->g.my_cellCentroid, ci);
    geo::copy(B, cc); geo::copy(B+3, xaxis); geo::copy(B+6, yaxis);

    ccvtx_lsinv_shape_function(order, t.c2ne[ci], xs, ys, B+9);
  }
}

void CellFromVtxsShapeFunction
::reinit (const Array3D<const Real>& p, const int lvl) {
  if (order < 2) return;
  const auto& t = mm->t;
  ompparfor for (int ci = 0; ci < ncell; ++ci) {
    // Approximate cell centroid. Only an approximation is needed b/c this just
    // sets the origin of a local coord system.
    const int ne = t.c2ne[ci];
    Vec cc = {0};
    for (int i = 0; i < ne; ++i) geo::axpy(1, p(lvl,t.get_c2v(ci,i)), cc);
    geo::scale(1.0/ne, cc);
    geo::normalize(cc);

    // Local x-y coord system.
    Vec xaxis, yaxis;
    Real xs[max_halo_slots], ys[max_halo_slots];
    map_vtxs_to_cell_plane(*mm, ci, xaxis, yaxis, xs, ys,
                           nullptr, scale_by_arc_length, cc, &p, lvl);

    Real* B = ho_coefs(ci);
    geo::copy(B, cc); geo::copy(B+3, xaxis); geo::copy(B+6, yaxis);

    ccvtx_lsinv_shape_function(order, t.c2ne[ci], xs, ys, B+9);
  }
}

namespace citest {

Real test_fn (const Int f, const Real x, const Real y, const Real z) {
  const Real r = std::sqrt(x*x + y*y + z*z);
  return ((1 + std::cos(f*M_PI*(x/r)))*
          square(1 + std::cos(f*M_PI*(y/r)))*
          (1 + std::sin(f*M_PI*(z/r))));
}

// Fill cell data. Comm is not needed because we fill all cells.
void fill_cc_data (const MpasMesh::ConstPtr& mm, const Int fac, const Ar2& f,
                   const int tq_order) {
  const auto fn = [&] (const Real x, const Real y, const Real z) {
    return test_fn(fac, x, y, z);
  };
  const auto ti = create_TestIntegrand(fn);
  const auto ncell = mm->t.ncell, nlvl = f.m();
  CellIntegrator cell_int(mm, tq_order);
  ompparfor for (Int ci = 0; ci < ncell; ++ci) {
    Real area, mass;
    cell_int.integrate(ci, ti, area, mass);
    for (Int li = 0; li < nlvl; ++li) f(li,ci) = mass/area;
  }
}

} // namespace citest

using namespace citest;

struct UnitTest {
  static constexpr Int tq_order = 8;
  const MpasMesh::ConstPtr mm;
  const siqp::sh::Mesh& m;
  const Comm::Ptr he;
  const Topo& t;
  const Vtxs& vtxs;
  const Int ncown;
  CellFromCellsModalShapeFunction ccsf;
  CellFromCellAndVtxsModalShapeFunction cvsf;
  CellFromVtxsShapeFunction cvnmsf;
  Ar1 q0_cc;
  Ar2 q_vtx;
  Ar2 shapes;

  UnitTest  (const MpasMesh::ConstPtr& mm_, const siqp::sh::Mesh& m_,
             const Comm::Ptr& he_)
    : mm(mm_), m(m_), he(he_), t(mm->t), vtxs(mm->v), ncown(mm->ncown)
  {}

  void fill_vtx_data (const Int fac) {
    ompparfor for (Int ci = 0; ci < ncown; ++ci) {
      const auto ne = t.c2ne[ci];
      for (Int i = 0; i < ne; ++i) {
        const auto vi = t.get_c2v(ci,i);
        Vec vtx; fill(vtxs, vi, vtx, true);
        q_vtx(i,ci) = test_fn(fac, vtx[0], vtx[1], vtx[2]);
      }
    }
  }

  // Compute shape coefs from cell-center data.
  void calc_shape_coefs_from_ccs () {
    if (ccsf.order == 1) return;
    Acr3 f_cc(1, 1, t.ncell, q0_cc.data());
    ompparfor for (Int ci = 0; ci < ncown; ++ci)
      ccsf.calc_shape(ci, f_cc, 0, 0, shapes(ci));
  }

  // Evaluate shape functions at cell vertices.
  void eval_shapes_at_vtxs () {
    Acr3 f_cc(1, 1, t.ncell, q0_cc.data());
    ompparfor for (Int ci = 0; ci < ncown; ++ci)
      ccsf.eval_shape_at_vtxs(ci, f_cc(0,0,ci), shapes(ci), 0, 0, q_vtx(ci));
  }

  bool
  test_shapes_against_coefs_from_vtxs (const CellFromVtxsShapeFunctionBase& sf,
                                       const int os, const Real tol) {
    if (sf.order == 1) return true;
    Acr4 f_vtx(max_vtx, 1, 1, ncown, q_vtx.data());
    bool ok = true;
    ompparfor for (Int ci = 0; ci < ncown; ++ci) {
      Real shape[sf::max_ho_ncoef];
      sf.calc_shape(ci, q0_cc(ci), f_vtx, 0, 0, shape);
      Real max_coef = 0;
      for (int i = 0; i < sf.ho_ncoef; ++i)
        max_coef = std::max(max_coef, std::abs(shape[i]));
      const auto* const shape_orig = shapes(ci);
      for (int i = 0; i < sf.ho_ncoef-os; ++i)
        if (std::abs(shape_orig[i] - shape[os+i]) > tol*max_coef)
          ok = false;
    }
    Real l = ok ? 0 : 1, g;
    he->all_reduce(Comm::MAX, 1, &l, &g);
    if (g > 0 && he->amroot())
      printf("test_shapes_against_coefs_from_vtxs FAILed\n");
    return ok;
  }

  // Compute shape coefs from vtx data.
  void calc_shape_coefs_from_vtxs (const CellFromVtxsShapeFunctionBase& sf) {
    if (sf.order == 1) return;
    Acr4 f_vtx(max_vtx, 1, 1, ncown, q_vtx.data());
    ompparfor for (Int ci = 0; ci < ncown; ++ci)
      sf.calc_shape(ci, q0_cc(ci), f_vtx, 0, 0, shapes(ci));
  }

  bool test_integration (const CellShapeFunction& sf) {
    Int nerr = 0;
    TriangleQuadrature tq;
    const auto& g = mm->g;
    ompparfor for (Int ci = 0; ci < ncown; ++ci) {
      // Since this is an Eulerian Voronoi polygon, integrating over it in one
      // pass would produce 0 for all but the first moment. To test integration
      // fully, break the polygon into triangles and accumulate over these.
      Real i1, i2;
      const auto ne = t.c2ne[ci];
      Real buf[9];
      Ar2 vo(3, 3, buf);
      geo::copy(vo(0), m.p(m.e(0,ci)));
      i1 = g.my_areaCell[ci]*q0_cc(ci);
      i2 = 0;
      Real basis[sf::max_halo_slots];
      for (Int tri = 1; tri < ne-1; ++tri) {
        geo::copy(vo(1), m.p(m.e(tri,ci)));
        geo::copy(vo(2), m.p(m.e(tri+1,ci)));
        if (sf.order > 1)
          i1 += sf.integrate_over_polygon(tq, vo, 3, shapes(ci), ci);
        sf.integrate_basis_over_polygon(tq, vo, 3, ci, basis);
        i2 += sf.integrate(q0_cc(ci), shapes(ci), basis);
      }
      i1 /= g.my_areaCell[ci];
      i2 /= g.my_areaCell[ci];
    }
    if (nerr) printf("test_integration %d\n", nerr);
    return nerr == 0;
  }

  void eval_sf_error (const std::string& name, const Int fac,
                      const CellShapeFunction& sf,
                      Real& l2_err, Real& mass_err) {
    const auto& g = mm->g;
    const auto f_ls_den = [&] (const Real x, const Real y, const Real z) -> Real {
      return square(test_fn(fac, x, y, z));
    };
    const auto ti_den = create_TestIntegrand(f_ls_den);
    CellIntegrator cell_int(mm, tq_order);
    Real l2_num = 0, l2_den = 0, mass_num = 0, mass_den = 0;
#   pragma omp parallel for reduction(+: l2_num, l2_den, mass_num, mass_den)
    for (Int ci = 0; ci < ncown; ++ci) {
      const auto f_reconstruct = [&] (const Real px, const Real py,
                                      const Real pz) -> Real {
        if (sf.order == 1) return q0_cc(ci);
        const Real* const cc = sf.ho_coefs(ci);
        const Real* const xaxis = cc + 3;
        const Real* const yaxis = cc + 6;
        const Vec p = {px, py, pz};
        Real x, y;
        calc_xy(cc, p, xaxis, yaxis, sf::scale_by_arc_length, x, y);
        const Real* mass_const = sf.mass_consts.empty() ? nullptr : sf.mass_consts(ci);
        return q0_cc(ci) + sf::eval_shape_fn(sf.order, shapes(ci),
                                             mass_const, x, y);
      };
      const auto f_ls_num = [&] (const Real x, const Real y, const Real z) -> Real {
        return square(f_reconstruct(x, y, z) - test_fn(fac, x, y, z));
      };
      const auto ti_num = create_TestIntegrand(f_ls_num);
      Real area, num, den;
      cell_int.integrate(ci, ti_num, area, num); l2_num += num;
      cell_int.integrate(ci, ti_den, area, den); l2_den += den;
      const auto ti_mass_num = create_TestIntegrand(f_reconstruct);
      cell_int.integrate(ci, ti_mass_num, area, num); mass_num += num;
      mass_den += q0_cc(ci)*g.my_areaCell[ci];
    }
    Real vl[] = {l2_num, l2_den, mass_num, mass_den}, vg[4];
    he->all_reduce(Comm::SUM, 4, vl, vg);
    l2_err = std::sqrt(vg[0]/vg[1]);
    mass_err = std::abs(vg[2] - vg[3])/vg[3];
  }

  bool run () {
    static const Real eps = std::numeric_limits<Real>::epsilon();
    static const Int facs[] = {4, 2, 1}, ntrial = Int(sizeof(facs)/sizeof(*facs));

    bool ok = true;

    q0_cc.reset(t.ncell);
    q_vtx.reset(max_vtx, ncown);
    shapes.reset(sf::max_ho_ncoef, ncown);

    for (const Int order : {1, 2, 3}) {
      ccsf.init(mm, m, he, order, tq_order, false, false);
      cvsf.init(mm, m, he, order, tq_order, false);
      // Pentagons don't permit a non-modal quadratic basis.
      cvnmsf.init(mm, he, order, tq_order, true);

      Real c_l2_err_prev = 0, v_l2_err_prev = 0, vnm_l2_err_prev = 0;
      for (Int trial = 0; trial < ntrial; ++trial) {
        if (order == 2 && trial == 1) {
          // Test reinit.
          cvnmsf.reinit(Array3D<const Real>(m.p.m(), 1, m.p.n(), m.p.data()), 0);
        }

        const auto fac = facs[trial];
        Real c_l2_err, v_l2_err, vnm_l2_err, mass_err;

        fill_cc_data(mm, fac, Ar2(1, t.ncell, q0_cc.data()), tq_order);
        calc_shape_coefs_from_ccs();
        eval_sf_error("ccsf", fac, ccsf, c_l2_err, mass_err);
        if (mass_err > 1e5*eps) { printf("ccsf mass err %1.3e\n", mass_err); ok = false; }
        if ( ! test_integration(ccsf)) ok = false;

        eval_shapes_at_vtxs();
        if ( ! test_shapes_against_coefs_from_vtxs(cvsf, 0, 1e-5)) ok = false;
        if (order == 2 && ! test_shapes_against_coefs_from_vtxs(cvnmsf, 1, 1e-2))
          ok = false;

        fill_vtx_data(fac);
        calc_shape_coefs_from_vtxs(cvsf);
        eval_sf_error("cvsf", fac, cvsf, v_l2_err, mass_err);
        if (mass_err > 1e5*eps) { printf("cvsf mass err %1.3e\n", mass_err); ok = false; }
        if ( ! test_integration(cvsf)) ok = false;

        calc_shape_coefs_from_vtxs(cvnmsf);
        eval_sf_error("cvnmsf", fac, cvnmsf, vnm_l2_err, mass_err);
        if (mass_err > 1e-3) { printf("cvnmsf mass err %1.3e\n", mass_err); ok = false; }
        if ( ! test_integration(cvnmsf)) ok = false;

        const auto calc_ooa = [&] (Real ed, Real gf) {
          return std::log(ed)/std::log(gf);
        };
        const Real gf = trial == 0 ? 0 : facs[trial-1]/facs[trial];
        const Real c_ooa = trial == 0 ? 0 : calc_ooa(c_l2_err_prev/c_l2_err, gf);
        const Real v_ooa = trial == 0 ? 0 : calc_ooa(v_l2_err_prev/v_l2_err, gf);
        const Real vnm_ooa = (trial == 0 ? 0 :
                              calc_ooa(vnm_l2_err_prev/vnm_l2_err, gf));
        if (trial == ntrial-1) {
          const Real tol = 0.95;
          const bool c_ok = c_ooa >= tol*order;
          const bool v_ok = v_ooa >= tol*order;
          if ( ! (c_ok && v_ok)) ok = false;
        }
        if (he->amroot()) {
          if (trial == 0 && order == 1)
            printf("ut_csf> order fac    ccsf l2 (ooa)         cvsf"
                   "             cvnmsf\n");
          printf("ut_csf>     %d  %2d %10.3e (%5.3f) %10.3e (%5.3f)",
                 order, fac, c_l2_err, c_ooa, v_l2_err, v_ooa);
          printf(" %10.3e (%5.3f)\n", vnm_l2_err, vnm_ooa);
        }

        c_l2_err_prev = c_l2_err;
        v_l2_err_prev = v_l2_err;
        vnm_l2_err_prev = vnm_l2_err;
      }
    }

    return ok;    
  }
};

bool unittest_CellShapeFunction (const MpasMesh::ConstPtr& mm,
                                 const siqp::sh::Mesh& m, const Comm::Ptr& he) {
  UnitTest t(mm, m, he);
  return t.run();
}
