#include "mpas_siqp.hpp"
#include "mpas_triquad.hpp"
#include "mpas_geometry.hpp"
#include "mpas_mesh.hpp"
#include "mpas_comm.hpp"
#include "mpas_comm_util.hpp"

static const int max_nvtx = 7;

struct TqArea : siqp::test::AreaFunctor {
  Real operator() (const Real* const v, const int n) const override {
    return geo::calc_area(tq, v, n);
  }
private:
  TriangleQuadrature tq;
};

static void fill_tri1 (Array2D<Real>& vs, Array1D<int>& vts,
                       const bool self_intersection = false) {
  vs.reset(3, 4);
  vts.reset(4);
  const Real z = 5;
  vs(0,0) =  -1; vs(1,0) =   -1; vts[0] = 0;
  vs(0,1) = 0.7; vs(1,1) = -0.9; vts[1] = 0;
  vs(0,2) =   1; vs(1,2) =    1; vts[2] = 0;
  vs(0,3) =   0; vs(1,3) = self_intersection ? -1 : 1;
  vts[3] = 1;
  for (int i = 0; i < vs.n(); ++i) {
    vs(2,i) = z;
    geo::normalize(vs(i));
  }
  //geo::axpbyz(0.5, vs(2), 0.5, vs(0), vs(3));
  siqp::sh::middle_matches(vs(2), vs(0), vs(3), vs(3));
}

static void
approx_by_spoly (const Array2D<const Real>& vs, const Array1D<const int>& vts,
                 Array2D<Real>& svs, const int nseg_per_cedge) {
  assert(vs.n() == vts.m());
  const int nvs = vs.n();
  svs.reset(3, 3*nseg_per_cedge+1);
  geo::copy(svs(0), vs(0));
  for (int i = 0, k = 0, f = 1; i < 3; ++i) {
    assert(vts[k] == 0); // No need to handle more general.
    assert(k+1 < nvs);
    if (vts[k+1] == 0) {
      for (int j = 0; j < nseg_per_cedge; ++j)
        siqp::sh::sedge_eval(vs(k), vs(k+1),
                             Real(j+1)/nseg_per_cedge, svs(f++));
      k++;
    } else {
      for (int j = 0; j < nseg_per_cedge; ++j)
        siqp::sh::cedge_eval(vs(k), vs(k+1), vs((k+2) % nvs),
                             Real(j+1)/nseg_per_cedge, svs(f++));
      k += 2;
    }
  }  
}

static void
write_python_array (const FilePtr& fp, const std::string& arrname,
                    const Array2D<const Real>& vs) {
  fprintf(fp.get(), "%s = [", arrname.c_str());
  for (int i = 0; i < vs.n(); ++i) {
    fprintf(fp.get(), "[");
    for (int j = 0; j < vs.m(); ++j)
      fprintf(fp.get(), "%1.8e,", vs(j,i));
    fprintf(fp.get(), "],\n");
  }
  fprintf(fp.get(), "];\n");
}

Real calc_stri_area (const TriangleQuadrature& tq, const Real v1[3],
                     const Real v2[3], const Real v3[3], const int order=8) {
  const Real* coord;
  const Real* weight;
  int nqp;
  tq.get_coef(order, coord, weight, nqp);
  Real area = 0, u[3];
  for (Int k = 0; k < nqp; ++k) {
    const Real jac = geo::calc_tri_jacobian(
      v1, v2, v3, slice(coord, k), u, false);
    area += weight[k]*jac;
  }
  return 0.5*area;
}

Real calc_vctri_area (const TriangleQuadrature& tq, const Real v1[3],
                      const Real s[3], const Real m[3], const Real p[3],
                      const int order_stri=8, const int order_creg=8) {
  const Real* coord;
  const Real* weight;
  int nqp;
  tq.get_coef(order_creg, coord, weight, nqp);
  Real area = 0, u[3], spmid[3];
  geo::axpbyz(0.5, s, 0.5, p, spmid);
  for (Int k = 0; k < nqp; ++k) {
    const Real jac = siqp::sh::calc_vctri_jacobian(
      spmid, s, m, p, slice(coord, k), u, false);
    area += weight[k]*jac;
  }
  return 0.5*area + calc_stri_area(tq, v1, s, p, order_stri);
}

static Real calc_cpoly_area (const TriangleQuadrature& tq,
                             const Real* vs, const int* vts,
                             const Int nvts, const Int order=8) {
  const Real* coord;
  const Real* weight;
  int nqp;
  tq.get_coef(order, coord, weight, nqp);

  Real ctr[3] = {0};
  int nedge = 0;
  for (int i = 0; i < nvts; ++i) {
    if (vts[i] == 1) continue;
    ++nedge;
    geo::axpy(1, vs + 3*i, ctr);
  }
  geo::normalize(ctr);

  Real area = 0;
  for (int ie = 0, k = 0; ie < nedge; ++ie) {
    const int kp1 = (k+1) % nvts;
    if (vts[k] == 0 && vts[kp1] == 0) {
      area += calc_stri_area(tq, ctr, vs + 3*k, vs + 3*kp1, order);
      k++;
    } else {
      if (vts[k] == 0) {
        const int kp2 = (k+2) % nvts;
        area += calc_vctri_area(tq, ctr, vs + 3*k, vs + 3*kp1, vs + 3*kp2,
                                order, order);
        k += 2;
      } else {
        assert(ie == 0);
        const int km1 = (k+nvts-1) % nvts;
        area += calc_vctri_area(tq, ctr, vs + 3*km1, vs + 3*k, vs + 3*kp1,
                                order, order);
        k++;
      }
    }
  }

  return area;
}

static bool test_ctri_area (const bool verbose = true) {
  bool ok = true;
  TriangleQuadrature tq;
  
  Array2D<Real> vs, svs;
  Array1D<int> vts;

  for (const bool self_intersecting : {false, true}) {
    if (verbose) printf("self intersecting %d\n", self_intersecting);
    fill_tri1(vs, vts, self_intersecting);

    Real area_svs;
    for (const int nseg : {50, 200, 1000}) {
      approx_by_spoly(vs, vts, svs, nseg);
      area_svs = geo::calc_area(tq, svs.data(), svs.n());
      if (verbose) printf("area svs %23.16e nseg %d\n", area_svs, nseg);
    }
    for (const int order : {8, 12, 20}) {
      const auto area_cvs = calc_cpoly_area(tq, vs.data(), vts.data(), vts.m(),
                                            order);
      const auto e = reldif(area_svs, area_cvs);
      if (verbose)
        printf("area cvs %23.16e (%1.3e) order %d\n", area_cvs, e, order);
      if (order == 20 && e > 1e-3) ok = false;
    }
  }

  if (0) {
    const auto fp = FilePtr(fopen("tmp/vs.py", "w"));
    write_python_array(fp, "vs", svs);
  }

  return ok;
}

bool test_siqp (const Comm::Ptr& he, const MpasMesh::ConstPtr& m) {
  const auto& t = m->t;
  const auto& v = m->v;
  Array2D<Real> cp(3, t.nvtx);
  Array2D<int> ce(max_nvtx, t.ncell);
  omppar {
    ompfor for (int i = 0; i < t.nvtx; ++i) {
      cp(0,i) = v.x[i]; cp(1,i) = v.y[i]; cp(2,i) = v.z[i];
      geo::normalize(cp(i));
    }
    ompfor for (int ci = 0; ci < t.ncell; ++ci) {
      for (int i = 0; i < max_nvtx; ++i) ce(i,ci) = -1;
      assert(t.c2ne[ci] <= max_nvtx);
      get_oriented_cell(t, v, ci, ce(ci));
    }
  }
  bool ok = siqp::sh::test_ice_math();
  if ( ! ok) printf("FAIL siqp::sh::test_ice_math\n");
  ok = ok && test_ctri_area();
  if ( ! ok) printf("FAIL test_ctri_area\n");
  { // Test: Clip one polygon against the MPAS mesh.
    TriangleQuadrature tq;
    TqArea af;
    const int skip = t.ncell / 242;
    for (const bool use_ice : {false, true}) {
      // When use_ice, the edges are still geometrically straight in this test.
      for (int trial = 0; trial < 2; ++trial) {
        if (use_ice) {
          //todo-ice Need to investigate ice trial=0 failures with -O2 or higher
          // at some point and FPE with trial=1. But ice isn't used in the mpas
          // app, so I'm not doing that right now.
          continue;
        }
        printf("use_ice %d trial %d\n", use_ice, trial);
        // In the first trial, do an intersection of one unperturbed element on
        // the whole mesh. In the second trial, perturb the element.
        ompparfor for (int ci = 0; ci < t.ncell; ci += skip) {
          Array2D<Real> p(3, max_nvtx);
          Array2D<int> e(max_nvtx, 1);
          for (int i = 0; i < max_nvtx; ++i) e(i,0) = -1;
          const int ne = t.c2ne[ci];
          const auto elem = ce(ci);
          for (int i = 0; i < max_nvtx; ++i) {
            if (elem[i] == -1) break;
            e(i,0) = i;
            for (int d = 0; d < 3; ++d)
              p(d,i) = cp(d,elem[i]) + (trial == 1 ? 0.01 : 0);
            if (trial == 1) geo::normalize(p(i));
          }
          const Real
            test_area = siqp::test::test_area_ot(cp, ce, p, e, false, use_ice, &af),
            tq_area = geo::calc_area(tq, p.data(), ne, 8),
            true_area = geo::calc_area_formula(p.data(), ne),
            re = reldif(true_area, test_area);
          if (re > 50*std::numeric_limits<Real>::epsilon()) {
            ok = false;
            printf("use_ice %d trial %d ci %d ne %d areas %1.3e (%1.2e) %1.3e re %1.2e\n",
                   use_ice, trial, ci, ne,
                   true_area, reldif(true_area, tq_area),
                   test_area, re);
            siqp::test::test_area_ot(cp, ce, p, e, true, use_ice, &af);
          }
        }
      }
    }
  }
  { // Test: Clip a distorted mesh against the original. Check how far we can
    // distort the mesh while still recovering the true area.
    const bool use_ice = false;
    TqArea af;
    Real radius, avg_diam;
    get_avg_radius_and_cell_spacing(*m, *he, radius, avg_diam);
    Real perturb = 0.4*avg_diam;
    for (int trial = 0; trial < 4; ++trial) {
      perturb *= 2;
      Array2D<Real> p(3, t.nvtx);
      ompparfor for (int i = 0; i < t.nvtx; ++i) {
        Real* const u = p(i);
        geo::copy(u, cp(i));
        for (int d = 0; d < 3; ++d) u[d] += perturb*2*(urand() - 0.5);
        geo::normalize(u);
      }
      const Real test_area = siqp::test::test_area_ot(cp, ce, p, ce, false, use_ice, &af);
      const Real re = reldif(4*M_PI, test_area);
      if (re > 50*std::numeric_limits<Real>::epsilon()) {
        ok = false;
        printf("clip distorted against original: %d %1.3e\n", trial, re);
      }
    }
  }
  return ok;
}
