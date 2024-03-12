#include "slmm_defs.hpp"
#include "slmm_mesh.hpp"
#include "slmm_io.hpp"
#include "slmm_nla.hpp"
#include "slmm_spf.hpp"
#include "slmm_time_int.hpp"
#include "slmm_gallery.hpp"
#include "slmm_debug.hpp"
#include "slmm_fit_extremum.hpp"
#include "slmm_gll.hpp"
#include "slmm_islet.hpp"
#include "slmm_basis_reduced.hpp"

using namespace slmm;

struct Command {
  enum Enum {
    test_make_cubedsphere, test_make_gll_mesh, test_make_gll_subcell_mesh,
    test_gll, test_gll_2d, test_time_int, test_qp_limiter, test_face_tree,
    test_spf, test_fit_extremum, test_nla, islet_compute, test_mass_matrix
  };
  static Enum from_string (const std::string& s) {
    if (s == "test_make_cubedsphere") return test_make_cubedsphere;
    if (s == "test_make_gll_mesh") return test_make_gll_mesh;
    if (s == "test_make_gll_subcell_mesh") return test_make_gll_subcell_mesh;
    if (s == "test_gll") return test_gll;
    if (s == "test_gll_2d") return test_gll_2d;
    if (s == "test_time_int") return test_time_int;
    if (s == "test_qp_limiter") return test_qp_limiter;
    if (s == "test_face_tree") return test_face_tree;
    if (s == "test_spf") return test_spf;
    if (s == "test_fit_extremum") return test_fit_extremum;
    if (s == "test_nla") return test_nla;
    if (s == "islet_compute") return islet_compute;
    if (s == "test_mass_matrix") return test_mass_matrix;
    throw std::runtime_error(s + " is not a command.");
  }
};

struct Input {
  Command::Enum command;
  Int ne, np;
  bool write_matlab, quiet;
  std::string fn_pre_out;

  Input(Int argc, char** argv);
  void print(std::ostream& os) const;
};

static Int test_make_cubedsphere (const Input& in) {
  AVec3s cp;
  AIdxs c2n;
  mesh::make_cubedsphere_mesh(cp, c2n, in.ne);
  Int nerr = 0;
  {
    const Int ne = mesh::impl::check_elem_normal_against_sphere(cp, c2n);
    if (ne) std::cerr << "FAIL: check_elem_normal_against_sphere\n";
    nerr += ne;
  }
  {
    AIdxs c2e, e2n;
    mesh::impl::make_c2e_from_c2n(in.np, c2n, c2e, e2n);
    Int ne = 0;
    // Every edge has two cells, and each cell is a quad.
    if (nslices(e2n) != 4/2*nslices(c2n)) {
      ++ne;
      std::cerr << "FAIL: make_c2e_from_c2n\n";
    }
    nerr += ne;
  }
  if (in.write_matlab)
    write_matlab("cm", cp, c2n);
  {
    bool success = true;
    for (Int ci = 0, ncell = 6*in.ne*in.ne; ci < ncell; ++ci) {
      const auto& cell = slice(c2n,ci);
      // Get a point in the cell, approximately its cell center, on the sphere.
      Real xyz[3] = {0};
      for (Int vi = 0; vi < 4; ++vi) {
        const auto& p = slice(cp, cell[vi]);
        for (Int j = 0; j < 3; ++j) xyz[j] += p[j];
      }
      for (Int j = 0; j < 3; ++j) xyz[j] /= 4;
      Real den = 0;
      for (Int j = 0; j < 3; ++j) den += siqk::square(xyz[j]);
      den = std::sqrt(den);
      for (Int j = 0; j < 3; ++j) xyz[j] /= den;
      // Check whether we recover the cell index.
      const Int ci_calc = mesh::get_cell_idx(in.ne, 0, nullptr,
                                             xyz[0], xyz[1], xyz[2]);
      if (ci_calc != ci) {
        ++nerr;
        success = false;
      }
    }
    if ( ! success) std::cerr << "FAIL: get_cell_idx\n";
  }
  return nerr;
}

static Int test_gll (const Input& in) {
  Int nerr = 0;
  const Real tol = 1e2*std::numeric_limits<Real>::epsilon();
  GLL gll;
  const Real* x, * wt;
  for (Int np = 2; np <= GLL::np_max; ++np) {
    gll.get_coef(np, x, wt);
    if ( ! x) continue;
    Real sum = 0;
    for (Int i = 0; i < np; ++i)
      sum += wt[i];
    if (std::abs(2 - sum) > tol) {
      std::cerr << "test_gll " << np << ": 2 - sum = " << 2 - sum << "\n";
      ++nerr;
    }
    for (Int j = 0; j < np; ++j) {
      Real gj[GLL::np_max];
      gll.eval(np, x[j], gj);
      for (Int i = 0; i < np; ++i) {
        if (j == i) continue;
        if (std::abs(gj[i]) > tol) {
          std::cerr << "test_gll " << np << ": gj[" << i << "] = " << gj[i] << "\n";
          ++nerr;
        }
      }
    }
    {
      // Integrate to check weights. Use smallest np suited to np=16.
      const Int qn = 10;
      const Real* qx, * qwt;
      gll.get_coef(qn, qx, qwt);
      Real sum[GLL::np_max] = {0};
      for (Int k = 0; k < qn; ++k) {
        Real gk[GLL::np_max];
        gll.eval(np, qx[k], gk);
        for (Int j = 0; j < np; ++j)
          sum[j] += gk[j]*qwt[k];
      }
      for (Int j = 0; j < np; ++j) {
        const Real analytical = wt[j];
        const Real rd = reldif(analytical, sum[j]);
        if (rd > (np <= 13 ? 15 : 50)*std::numeric_limits<Real>::epsilon()) {
          printf("%d %d %1.3e %1.3e %1.1e\n", np, j, analytical, sum[j], rd);
          ++nerr;
        }
      }
    }
  }
  for (Int np = 2; np <= GLL::np_max; ++np) {
    gll.get_coef(np, x, wt);
    if ( ! x) continue;
    Real a[] = {-0.9, -0.7, -0.3, 0.1, 0.2, 0.4, 0.6, 0.8};
    const Real delta = std::sqrt(std::numeric_limits<Real>::epsilon());
    for (size_t ia = 0; ia < sizeof(a)/sizeof(Real); ++ia) {
      Real gj[GLL::np_max], gjp[GLL::np_max], gjm[GLL::np_max];
      gll.eval_derivative(np, a[ia], gj);
      gll.eval(np, a[ia] + delta, gjp);
      gll.eval(np, a[ia] - delta, gjm);
      for (Int i = 0; i < np; ++i) {
        const Real fd = (gjp[i] - gjm[i])/(2*delta);
        const Real d = std::abs(fd - gj[i]);
        Real tol = delta*std::abs(gjp[i]);
        if (np >= 5) tol *= 5;
        if (d >= tol) {
          printf("%d %d %1.3e %1.3e %1.1e %1.1e\n",
                 np, i, fd, gj[i], d, tol);
          ++nerr;
        }
      }
    }
  }
  return nerr;
}

static Int test_gll_2d (const Input& in) {
  Int nerr = 0;
  static const Real coords[4][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
  static const Int tris[2][3] = {{0, 1, 2}, {0, 2, 3}};

  // See comment below re: tq_order.
  for (int np = 2; np <= std::min<Int>(11, GLL::np_max); ++np) {
    GLL gll;
    const Real* gll_x, * gll_wt;
    gll.get_coef(np, gll_x, gll_wt);

    siqk::TriangleQuadrature tq;
    siqk::RawConstVec3s tq_bary;
    siqk::RawConstArray tq_w;
    // This is enough for np <= 11, but it's not enough for np = 12.
    const int tq_order = 20;
    tq.get_coef(tq_order, tq_bary, tq_w);

    Real numerical[GLL::np_max*GLL::np_max] = {0};
    for (int trii = 0; trii < 2; ++trii) {
      const auto tri = tris[trii];
      for (Int q = 0; q < len(tq_w); ++q) {
        Real coord[2];
        siqk::PlaneGeometry::bary2coord(
          coords[tri[0]], coords[tri[1]], coords[tri[2]],
          slice(tq_bary, q),
          coord);
        Real gi[GLL::np_max], gj[GLL::np_max];
        gll.eval(np, coord[0], gi); 
        gll.eval(np, coord[1], gj);

        for (int i = 0; i < np; ++i)
          for (int j = 0; j < np; ++j)
            numerical[i*np + j] += 0.5 * tq_w[q] * gi[i] * gj[j];
      }
    }

    Real tol = 5*std::numeric_limits<Real>::epsilon();
    if (np >= 5) tol *= 5;
    if (np >= 7) tol *= 10;
    for (int i = 0; i < np; ++i)
      for (int j = 0; j < np; ++j) {
        const Real analytical = 0.25 * gll_wt[i] * gll_wt[j];
        const Real rd = reldif(analytical, numerical[np*i + j]);
        if (rd > tol) {
          printf("%d %d %d %1.3e %1.3e %1.1e\n",
                 np, i, j,
                 analytical, numerical[np*i + j], rd);
          ++nerr;
        }
      }
  }

  return nerr;
}

static Int test_make_gll_mesh (const Input& in,
                               const bool test_make_gll_subcell_mesh=false) {
  AVec3s geo_p, cgll_p;
  AIdxs geo_c2n, cgll_c2n, cgll_io_c2n;
  const Int np = test_make_gll_subcell_mesh ? 2 : in.np;
  if (test_make_gll_subcell_mesh)
    mesh::make_cubedsphere_subcell_mesh(in.ne, in.np, GLL(), geo_p, geo_c2n);
  else
    mesh::make_cubedsphere_mesh(geo_p, geo_c2n, in.ne);
  Int nerr = 0;
  mesh::make_cgll_from_geo(geo_p, geo_c2n, np, GLL(), cgll_p, cgll_c2n);
  mesh::make_io_cgll_from_internal_cgll(cgll_p, cgll_c2n, cgll_io_c2n);
  { // Clip the mesh against itself and get the total area.
    const Real
      area = siqk::test::test_area_ot<siqk::SphereGeometry>(
        a2Vec3s(cgll_p), a2Idxs(cgll_io_c2n), a2Vec3s(cgll_p), a2Idxs(cgll_io_c2n)),
      true_area = 4*M_PI,
      re = std::abs(area - true_area)/true_area;
    if (re >= 1e-10) {
      fprintf(stderr, "true area %1.4e mesh area %1.4e relerr %1.4e\n",
              true_area, area, re);
      ++nerr;
    }
  }
  {
    const Int ne = mesh::impl::check_elem_normal_against_sphere(
      cgll_p, cgll_io_c2n);
    if (ne) std::cerr << "FAIL: check_elem_normal_against_sphere\n";
    nerr += ne;
  }
  {
    AIdxArray dglln2cglln;
    AIdxs dgll_c2n;
    mesh::make_dgll_from_cgll(cgll_p, cgll_c2n, dglln2cglln, dgll_c2n);
    const Int np2 = szslice(cgll_c2n);
    Int ne = 0;
    for (Int ci = 0; ci < nslices(cgll_c2n); ++ci) {
      const auto cgll_cell = slice(cgll_c2n, ci);
      for (Int ni = 0; ni < siqk::square(np); ++ni)
        if (dglln2cglln[ci*np2 + ni] != cgll_cell[ni])
          ++ne;
    }
    if (ne) {
      nerr += ne;
      std::cerr << "FAIL: make_dgll_from_cgll\n";
    }
  }
  if ( ! in.fn_pre_out.empty()) {
    io::NetcdfWriter ncw(cgll_p, cgll_io_c2n, in.fn_pre_out + ".g");
    ncw.add_nodal_field("x");
    ncw.end_definition();
    const Int n = nslices(cgll_p);
    std::vector<Real> x(n), lat(n), lon(n);
    for (Int i = 0; i < n; ++i) {
      const auto p = slice(cgll_p, i);
      xyz2ll(p[0], p[1], p[2], lat[i], lon[i]);
    }
    ncw.advance_time_to(1);
    gallery::InitialCondition::init(
      gallery::InitialCondition::CosineBells,
      nslices(cgll_p), lat.data(), lon.data(), x.data());
    ncw.write_field("x", x.data());
    ncw.advance_time_to(1.5);
    gallery::InitialCondition::init(
      gallery::InitialCondition::SlottedCylinders,
      nslices(cgll_p), lat.data(), lon.data(), x.data());
    ncw.write_field("x", x.data());
    ncw.advance_time_to(2.5);
    gallery::InitialCondition::init(
      gallery::InitialCondition::CorrelatedCosineBells,
      nslices(cgll_p), lat.data(), lon.data(), x.data());
    ncw.write_field("x", x.data());
  }
  if (in.write_matlab) {
    write_matlab("cm", geo_p, geo_c2n);
    write_matlab("m", cgll_p, cgll_io_c2n);
    write_matlab("gll", cgll_p, cgll_c2n);
  }
  { // Test get_adjacent_cells.
    AIdxArray geo_c2cnbrs_ptr, geo_c2cnbrs;
    mesh::get_adjacent_cells(geo_c2n, geo_c2cnbrs_ptr, geo_c2cnbrs);
    for (Int i = 0; i < nslices(geo_c2cnbrs_ptr) - 1; ++i) {
      for (Int j = geo_c2cnbrs_ptr[i]; j < geo_c2cnbrs_ptr[i+1]; ++j) {
        // Check that cell i and cell geo_c2cnbrs[j] have >= 1 node in common.
        bool found = false;
        for (Int k0 = 0; k0 < 4; ++k0) {
          for (Int k1 = 0; k1 < 4; ++k1)
            if (geo_c2n(geo_c2cnbrs[j],k0) == geo_c2n(i,k1)) {
              found = true;
              break;
            }
          if (found) break;
        }
        if ( ! found) ++nerr;
      }
    }
  }
  return nerr;
}

static Int test_time_int (const Input& in) {
  return timeint::test::test_ark( ! in.quiet);
}

static Int test_qp_limiter (const Input& in) {
  const Int n = square(in.np);
  GLL gll;
  const Real* gx, * gw;
  gll.get_coef(in.np, gx, gw);

  auto test_fn = [=] (const Real x, const Real y) -> Real {
    return std::sin(0.3*2*M_PI*(x + 0.5))*std::cos(0.4*2*M_PI*(y + 0.4));
  };
  const Real f_regression[] = {
    0.0523002, 0.109336, 0.164448, 0.122502, -0.700229, 0.201736, 0.5, 0.409943,
    0.5, 0.0504606, -0.41464, -0.0606496, 0.5, 0.0107386, -0.805338, -0.184218 };

  std::vector<Real> w(n), flo(n, -0.9), fhi(n, 0.5), f(n);
  for (Int i = 0, k = 0; i < in.np; ++i)
    for (Int j = 0; j < in.np; ++j, ++k) {
      w[k] = gw[i]*gw[j];
      f[k] = test_fn(gx[j], gx[i]);
    }
  Real b = 0;
  for (Int i = 0; i < n; ++i)
    b += w[i]*f[i];
  Int nerr = 0;
  for (int trial = 0; trial < 3; ++trial) {
    std::vector<Real> f_lim(f);
    if (trial == 2) { flo[0] = -100; fhi[0] = 100; }
    const Int info = spf::solve_1eq_bc_qp(n, w.data(), w.data(), b, flo.data(),
                                          fhi.data(), trial >= 1, f.data(),
                                          f_lim.data(), 3);
    if (trial <= 1) {
      if (info != 1) ++nerr;
      Real err = 0;
      for (Int i = 0; i < n; ++i) err += square(f_lim[i] - f_regression[i]);
      err = std::sqrt(err);
      if (err > 1e-6) ++nerr;
    } else {
      if (info != 0) ++nerr;
    }
  }
  return nerr;
}

static Int test_face_tree (const Input& in) {
  static const Int nes[] = {1, 2, 3, 4, 5, 10, 11, 20, 33, 120};
  static const Int nss[] = {0, 1, 0, 2, 0, 3,  0,  1,  2,  1};
  Int nerr = 0;
  for (size_t i = 0; i < sizeof(nes)/sizeof(Int); ++i) {
    AIdxArray tree;
    if ( ! mesh::make_cubedsphere_tree_over_cells(nes[i], tree, nss[i], true))
      ++nerr;
  }
  return nerr;
}

static Int test_spf (const Input& in) {
  Int nerr = 0;
  static const Int nes[] = {1, 2, 3, 4, 5, 10, 11, 20, 33, 120};
  for (size_t i = 0; i < sizeof(nes)/sizeof(Int); ++i)
    nerr += spf::test(nes[i]);
  nerr += spf::test();
  return nerr;
}

static Int test_nla (const Input& in) {
  BlockMatrix<>::test();
  return test_solve_kkt() + test_form_ls_op() + QrFac::unittest(11);
}

static Real calc_gll_area (const Basis& b, const Int np, const Vec3s& p, const Int nv) {
  static const Int e[] = {0,1,2,3};
  const Real* xnode, *wt;
  b.get_x(np, xnode);
  b.get_w(np, wt);
  Real area = 0;
  for (Int i = 0; i < np; ++i) {
    for (Int j = 0; j < np; ++j) {
      Real J[6], tmp[3];
      siqk::sqr::impl::calc_Jacobian(p, e, xnode[i], xnode[j], J);
      siqk::SphereGeometry::cross(J, J+3, tmp);
      const Real jac = std::sqrt(siqk::SphereGeometry::norm2(tmp));
      area += wt[i]*wt[j]*jac;
    }
  }
  return area/4;
}

static Int test_fit_extremum (const Input& in) {
  Int nerr = 0;
  Real y_gll[GLL::np_max*GLL::np_max], qec[9];
  const Real tol = 1e3*std::numeric_limits<Real>::epsilon();

  GLL gll;
  const Real* gx, * gw;
  gll.get_coef(in.np, gx, gw);

  FitExtremum qe(in.np);
  Real min, max;
  bool use;

  const Real fc[9] = {0, 0, 0, 3, -1.5, 4, -2, 0.5, 1};
  const auto f = [&] (Real x, Real y) {
    return fc[3]*x*x + fc[4]*x*y + fc[5]*y*y + fc[6]*x + fc[7]*y + fc[8];
  };
  for (Int i = 0; i < in.np; ++i)
    for (Int j = 0; j < in.np; ++j)
      y_gll[in.np*i + j] = f(gx[j], gx[i]);
  qe.calc(y_gll, min, max, use, qec);
  siqk::prarr("fc", fc, 9);
  siqk::prarr("qec", qec, 9);
  for (Int i = 0; i < 9; ++i)
    if (std::abs(fc[i] - qec[i]) > tol)
      ++nerr;

  for (int trial = 0; trial < 2; ++trial) {
    Real dx, dy;
    switch (trial) {
    case 0: dx = -0.8; dy = 0.4; break;
    case 1: dx = 1.1; dy = 0.9; break;
    }
    const auto g = [&] (Real x, Real y) {
      x += dx;
      y += dy;
      return 2 - (2*x*x - 1*x*y + 3*y*y);
    };
    for (Int i = 0; i < in.np; ++i)
      for (Int j = 0; j < in.np; ++j)
        y_gll[in.np*i + j] = g(gx[j], gx[i]);
    qe.calc(y_gll, min, max, use, qec);
    siqk::prarr("qec", qec, 9);
    pr(puf(use) pu(min) pu(max-2));
    if (( ! use && trial == 0) || (use && trial == 1)) ++nerr;
    if (use && std::abs(max - 2) > tol) ++nerr;
  }

  return nerr;
}

Int islet_compute (const Input& in) {
  islet::unittest_Nodes();
  GLL gll;
  const Int nerr = Basis::compute_and_print_weights(gll, false, true);
  printf("islet::GllOffsetNodal:\n");
  islet::GllOffsetNodal igon;
  Basis::compute_and_print_weights(igon);
  printf("islet::GllNodal:\n");
  islet::GllNodal ign;
  Basis::compute_and_print_weights(ign);
  printf("islet::FreeNodal:\n");
  islet::FreeNodal ifn;
  Basis::compute_and_print_weights(ifn);
  printf("UniformNodeReduced:\n");
  UniformNodeReduced gr;
  Basis::compute_and_print_weights(gr);
  printf("UniformOffsetNodal:\n");
  islet::UniformOffsetNodal uon;
  Basis::compute_and_print_weights(uon);
  return nerr;
}

Int test_mass_matrix (const Input& in) {
  constexpr Real eps = std::numeric_limits<Real>::epsilon();
  Int nerr = 0;
  std::vector<Real> M(square(square((Int)GLL::np_max)));
  for (Int trial = 0; trial < 3; ++trial) {
    Basis::ConstPtr b;
    if (trial == 0) b = std::make_shared<GLL>();
    else if (trial == 1) b = std::make_shared<islet::GllNodal>();
    else b = std::make_shared<UniformNodeReduced>();
    for (Int np = 4; np <= 16; ++np) {
      if (np > 13 && np < 16) continue;
      if (np > GLL::np_max) continue;
      const Real tolf = np <= 8 ? 40 : np <= 10 ? 80 : 100;
      const Int np2 = square(np);
      for (Int nf = 1; nf <= std::min(8, np); ++nf) {
        const bool ok = Basis::compute_integrals_over_subcells_2d(
          *b, np, nf, M.data(), true);
        if ( ! ok) ++nerr;
        // Test that sum of each column equals w[i]*w[j].
        const Real* w;
        b->get_w(np, w);
        for (Int i = 0; i < np; ++i)
          for (Int j = 0; j < np; ++j) {
            Real sum = 0;
            for (Int k = 0; k < nf*nf; ++k)
              sum += M[np2*k + (i*np + j)];
            if (reldif(w[i]*w[j], sum) > tolf*eps) {
              pr(puf(trial) pu(np) pu(nf) pu(i) pu(j) pu(w[i]*w[j]) pu(sum)
                 pu(reldif(w[i]*w[j], sum)));
              ++nerr;
            }
          }
      }
    }
  }
  for (Int trial = 0; trial < 2; ++trial) {
    // Expensive, so be a little selective.
    islet::GllNodal br;
    Basis::ConstPtr bc;
    if (trial == 0) bc = std::make_shared<GLL>();
    else if (trial == 1) bc = std::make_shared<UniformNodeReduced>();
    for (const Int npr : {4, 6, 9, 11, 13}) {
      if (npr > 13 && npr < 16) continue;
      const Int npr2 = square(npr);
      for (const Int npc : {4, 8, 11, 13}) {
        if (npc > 13 && npc < 16) continue;
        const Int npc2 = square(npc);
        const bool ok = Basis::compute_mass_matrix_2d(
          br, npr, *bc, npc, M.data(), true);
        if ( ! ok) ++nerr;
        // Test sum of rows and sum of columns.
        const Real tolf = npr <= 10 && npc <= 10 ? 150 : 1000;
        const Real* w;
        br.get_w(npr, w);
        for (Int i = 0; i < npr; ++i)
          for (Int j = 0; j < npr; ++j) {
            const Real* Mrow = &M[npc2*(npr*i + j)];
            Real sum = 0;
            for (Int k = 0; k < npc2; ++k) sum += Mrow[k];
            if (reldif(w[i]*w[j], sum) > tolf*eps) {
              pr("row" pu(trial) pu(npr) pu(npc) pu(i) pu(j) pu(w[i]*w[j]) pu(sum)
                 pu(reldif(w[i]*w[j], sum)));
              ++nerr;
            }
          }
        bc->get_w(npc, w);
        for (Int i = 0; i < npc; ++i)
          for (Int j = 0; j < npc; ++j) {
            const Int col = npc*i + j;
            Real sum = 0;
            for (Int k = 0; k < npr2; ++k) sum += M[k*npc2 + col];
            if (reldif(w[i]*w[j], sum) > tolf*eps) {
              pr("col" pu(trial) pu(npr) pu(npc) pu(i) pu(j) pu(w[i]*w[j]) pu(sum)
                 pu(reldif(w[i]*w[j], sum)));
              ++nerr;
            }
          }
      }
    }
  }
  return nerr;
}

Input::Input (Int argc, char** argv)
  : command(Command::test_make_cubedsphere), ne(10), np(4),
    write_matlab(false), quiet(false)
{
  for (Int i = 1; i < argc; ++i) {
    const std::string& token = argv[i];
    if (eq(token, "-c", "--command")) command = Command::from_string(argv[++i]);
    else if (eq(token, "-ne")) ne = atoi(argv[++i]);
    else if (eq(token, "-np")) np = atoi(argv[++i]);
    else if (eq(token, "-q", "--quiet")) quiet = true;
    else if (eq(token, "-m", "--write-matlab")) write_matlab = true;
    else if (eq(token, "-o", "--output-prefix")) fn_pre_out = argv[++i];
  }

  if ( ! quiet) print(std::cout);
}

void Input::print (std::ostream& os) const {
  os << "command (-c) " << command << "\n"
     << "ne (-ne): " << ne << "\n"
     << "np (-np): " << np << "\n"
     << "write matlab (-m): " << write_matlab << "\n";
}

int main (int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  {
    Input in(argc, argv);
    Int nerr = 0;
    switch (in.command) {
    case Command::test_make_cubedsphere: nerr = test_make_cubedsphere(in); break;
    case Command::test_make_gll_mesh: nerr = test_make_gll_mesh(in); break;
    case Command::test_make_gll_subcell_mesh: nerr = test_make_gll_mesh(in, true); break;
    case Command::test_gll: nerr = test_gll(in); break;
    case Command::test_gll_2d: nerr = test_gll_2d(in); break;
    case Command::test_time_int: nerr = test_time_int(in); break;
    case Command::test_qp_limiter: nerr = test_qp_limiter(in); break;
    case Command::test_face_tree: nerr = test_face_tree(in); break;
    case Command::test_spf: nerr = test_spf(in); break;
    case Command::test_nla: nerr = test_nla(in); break;
    case Command::test_fit_extremum: nerr = test_fit_extremum(in); break;
    case Command::islet_compute: nerr = islet_compute(in); break;
    case Command::test_mass_matrix: nerr = test_mass_matrix(in); break;
    }
    std::cerr << (nerr ? "FAIL" : "PASS") << "ED\n";
  }
  Kokkos::finalize_all();
}
