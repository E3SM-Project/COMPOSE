#include <omp.h>
#include "siqk_intersect.hpp"
#include "mexutil.hpp"
using namespace siqk;

static void make_elems (const mexutil::ConstDenseMexMat& me, Idxs& e) {
  for (size_t i = 0; i < me.n; ++i)
    for (size_t j = 0; j < me.m; ++j)
      e(i,j) = static_cast<int>(me.a[me.m*i + j]) - 1;
}

static void merror (const std::string& msg) {
  Kokkos::finalize();
  mexErrMsgTxt(msg.c_str());
}

void mexFunction (int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {
  omp_set_num_threads(4);
  Kokkos::initialize();
  using namespace mexutil;
  std::string cmd = init_mex(nrhs, prhs);
  try {
    typedef PlaneGeometry geo;
    if (cmd == "inside") {
      if (nlhs != 1 || nrhs != 2) merror("in = inside(edge, points)");
      ConstDenseMexMat edge(prhs[0]);
      reqorexit(edge.m == 3 && edge.n == 2);
      ConstDenseMexMat points(prhs[1]);
      reqorexit(points.m == 3);
      DenseMexMat in(1, points.n);
      plhs[0] = in.ma;
      for (size_t i = 0; i < points.n; ++i) {
        double en[3];
        geo::edge_normal(edge.a, edge.a + 3, en);
        in.a[i] = geo::inside(points.a + points.m*i, edge.a,
                              const_cast<const double*>(en));
      }
    } else if (cmd == "intersect") {
      // Assumption: Intersection exists.
      if (nlhs != 1 || nrhs != 2)
        merror("points = intersect(edge, edges)");
      ConstDenseMexMat edge(prhs[0]);
      reqorexit(edge.m == 3 && edge.n == 2);
      ConstDenseMexMat edges(prhs[1]);
      DenseMexMat points(edge.m, edges.n), exists(1, edges.n);
      plhs[0] = points.ma;
      for (size_t i = 0; i < edges.n; ++i) {
        double en[3];
        geo::edge_normal(edge.a, edge.a + 3, en);
        geo::intersect(edges.a + 6*i, edges.a + 6*i + 3, edge.a,
                       const_cast<const double*>(en), points.a + points.m*i);
      }
    } else if (cmd == "clip_against_edge") {
      if (nlhs != 1 || nrhs != 2)
        merror("vo = clip_against_edge(edge, vi)");
      ConstDenseMexMat edge(prhs[0]);
      reqorexit(edge.m == 3 && edge.n == 2);
      ConstDenseMexMat vi(prhs[1]);
      reqorexit(vi.m == 3);
      Vec3s vo("vo", test::max_nvert, 3);
      int no;
      double en[3];
      geo::edge_normal(edge.a, edge.a + 3, en);
      sh::clip_against_edge<geo>(RawConstVec3s(vi.a, vi.n, vi.m), vi.n,
                                 vo, no, edge.a, const_cast<const double*>(en));
      DenseMexMat vom(vi.n, no);
      memcpy(vom.a, vo.ptr_on_device(), vi.n*no*sizeof(double));
      plhs[0] = vom.ma;
    } else if (cmd == "clip_against_poly") {
      if (nlhs != 1 || nrhs != 2)
        merror("vo = clip_against_poly(clip_polygon, vi)");
      ConstDenseMexMat mcp(prhs[0]);
      reqorexit(mcp.m == 3);
      ConstDenseMexMat vi(prhs[1]);
      reqorexit(vi.m == 3);
      RawConstVec3s cp(mcp.a, mcp.n, mcp.m);
      Vec3s cens("cens", nslices(cp), 3);
      for (int i = 0; i < nslices(cp); ++i)
        geo::edge_normal(slice(cp,i), slice(cp, (i + 1) % nslices(cp)), slice(cens,i));
      Vec3s vo("vo", test::max_nvert, 3), wrk("wrk", test::max_nvert, 3);
      int no;
      sh::clip_against_poly<geo>(cp, cens,
                                 RawConstVec3s(vi.a, vi.n, vi.m), vi.n,
                                 vo, no, wrk);
      DenseMexMat vom(vi.m, no);
      memcpy(vom.a, vo.ptr_on_device(), vi.m*no*sizeof(double));
      plhs[0] = vom.ma;
    } else if (cmd == "clip_against_poly_sphere") {
      if (nlhs != 1 || nrhs != 2)
        merror("vo = clip_against_poly_sphere(clip_polygon, vi)");
      ConstDenseMexMat mcp(prhs[0]);
      reqorexit(mcp.m == 3);
      ConstDenseMexMat vi(prhs[1]);
      reqorexit(vi.m == 3);
      RawConstVec3s cp(mcp.a, mcp.n, mcp.m);
      Vec3s cens("cens", nslices(cp), 3);
      for (int i = 0; i < nslices(cp); ++i)
        SphereGeometry::edge_normal(slice(cp,i), slice(cp, (i + 1) % nslices(cp)),
                                    slice(cens,i));
      Vec3s vo("vo", test::max_nvert, 3), wrk("wrk", test::max_nvert, 3);
      int no;
      sh::clip_against_poly<SphereGeometry>(cp, cens,
                                            RawConstVec3s(vi.a, vi.n, vi.m), vi.n,
                                            vo, no, wrk);
      DenseMexMat vom(vi.m, no);
      memcpy(vom.a, vo.ptr_on_device(), vi.m*no*sizeof(double));
      plhs[0] = vom.ma;
#if 0
    } else if (cmd == "test_area_ot") {
      // Test using oct-tree.
      if (nlhs != 1 || nrhs != 4)
        merror("area = test_area_ot(cp, ce, p, e)");
      ConstDenseMexMat mcp(prhs[0]);
      reqorexit(mcp.m == 3);
      ConstDenseMexMat mce(prhs[1]);
      ConstDenseMexMat mp(prhs[2]);
      reqorexit(mp.m == 3);
      ConstDenseMexMat me(prhs[3]);
      Array2D<const double> cp(3, mcp.n, mcp.a);
      Array2D<const double> p(3, mp.n, mp.a);
      Array2D<int> ce(mce.m, mce.n), e(me.m, me.n);
      make_elems(mce, ce);
      make_elems(me, e);
      DenseMexMat area(1, 1);
      plhs[0] = area.ma;
      area.a[0] = test::test_area_ot<geo>(cp, ce, p, e);
    } else if (cmd == "test_area_ot_sphere") {
      // Test using oct-tree.
      if (nlhs != 1 || nrhs != 4)
        merror("area = test_area_ot(cp, ce, p, e)");
      ConstDenseMexMat mcp(prhs[0]);
      reqorexit(mcp.m == 3);
      ConstDenseMexMat mce(prhs[1]);
      ConstDenseMexMat mp(prhs[2]);
      reqorexit(mp.m == 3);
      ConstDenseMexMat me(prhs[3]);
      Array2D<const double> cp(3, mcp.n, mcp.a);
      Array2D<const double> p(3, mp.n, mp.a);
      Array2D<int> ce(mce.m, mce.n), e(me.m, me.n);
      make_elems(mce, ce);
      make_elems(me, e);
      DenseMexMat area(1, 1);
      plhs[0] = area.ma;
      area.a[0] = test::test_area_ot<SphereGeometry>(cp, ce, p, e);
#endif
    } else {
      merror((string("Invalid function: ") + cmd).c_str());
    }
  } catch (const std::exception& e) {
    merror(e.what());
  }
  Kokkos::finalize();
}
