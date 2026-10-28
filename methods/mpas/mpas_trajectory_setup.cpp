#include "mpas_trajectory_impl.hpp"

#include <set>

void Trajectory::Impl
::setup_advect_vtxs () {
  if (Method::is_exact(method)) {
    if (method == Method::exact) {
      initial_step.reset(mm->t.nvtx);
      initial_step.set(1e-3);
    }
    return;
  }

  const auto& t = mm->t;
  const auto& dist = he->cdhalo();

  // Cell haloes.
  require(nhalo > 0 && nhalo <= max_nhalo);
  c2cslotsz = max_halo_slots[nhalo-1];
  c2c.reset(c2cslotsz, t.ncell);
  c2ccnt.reset(t.ncell);
  ompparfor for (int ci = 0; ci < t.ncell; ++ci) {
    if (he->parallel() && dist(ci) && (Int) dist(ci) <= 1) {
      c2ccnt[ci] = -1; // for use in asserts
      continue;
    }
    c2ccnt[ci] = augment_stencil(*mm, ci, 0, c2c(ci), nhalo, he->parallel());
    throw_if(c2ccnt[ci] > c2cslotsz, "c2ccnt is " << c2ccnt[ci]
             << " but limit is " << c2cslotsz << " for nhalo " << nhalo);
  }
  cc_min_nbr_dist2.reset(t.ncell); {
    const auto& ccs = mm->g.my_cellCentroid;
    ompparfor for (int ci = 0; ci < t.ncell; ++ci) {
      Real min_dist2 = 1e30;
      const auto ne = t.c2ne[ci];
      for (int j = 0; j < ne; ++j) {
        const auto cj = t.get_c2c(ci,j);
        if (cj == -2) continue;
        assert(cj >= 0 && cj < t.ncell);
        Vec p;
        geo::axpbyz(1, ccs(ci), -1, ccs(cj), p);
        min_dist2 = std::min(min_dist2, geo::norm2(p));
      }
      cc_min_nbr_dist2(ci) = min_dist2;
    }
  }

  // 3D edge normals.
  std::vector<Real> enmls(3*t.nedge);
  ecs.reset(3, t.nedge);
  angle_edge.reset(t.nedge);
  ompparfor for (int ci = 0; ci < t.ncell; ++ci) {
    const int ne = t.c2ne[ci];
    Vec cc; fill(mm->c, ci, cc, true);
    for (int i = 0; i < ne; ++i) {
      const int ei = t.get_c2e(ci,i);
      int sign = 1;
      if (mm->g.my_sign[ci*mm->t.max_vtx + i] < 0) {
        if (he->parallel() && dist(ci) && (Int) dist(ci) <= 1) {
          // This cell is on the boundary of my shared cells.
          const int co = get_other_cell(t, ci, ei);
          // If the other cell is owned or shared, this normal will be computed
          // by that cell.
          if (co >= 0) continue;
          // Otherwise, this edge, and not just this cell, is on the
          // boundary. Compute its normal and account for the sign.
          sign = -1;
        } else {
          // Will be calculated when ci is the other cell.
          continue; 
        }
      }
      Vec ec; fill(mm->e, ei, ec, true);
      geo::copy(ecs(ei), ec);
      // Use this to mimic MPAS's procedure. Sign check above means our 3D edge
      // normal vector is consistent with sign and angleEdge.
      angle_edge[ei] = mm->g.my_angleEdge[ei];
      // Outward edge normal.
      Real* const enml = &enmls[3*ei];
      geo::axpbyz(1, ec, -1, cc, enml);
      // Make it tangent to the sphere.
      geo::axpy(-geo::dot(enml, ec), ec, enml);
      geo::normalize(enml);
      if (sign < 0) geo::scale(sign, enml);
#ifndef NDEBUG
      // Edge tangent, though we're not using it.
      Real etan[3];
      geo::cross(enml, ec, etan);
      geo::normalize(etan);
      const Real tol = 100*std::numeric_limits<Real>::epsilon();
      assert(std::abs(geo::dot(enml,ec)) < tol);
      assert(std::abs(geo::dot(enml,etan)) < tol);
      assert(std::abs(geo::dot(etan,ec)) < tol);
#endif
    }
  }
#ifndef NDEBUG
  // Every edge normal should be filled.
  ompparfor for (int ei = 0; ei < t.nedge; ++ei) {
    Real* const enml = &enmls[3*ei];
    const Real nrm = std::sqrt(geo::norm2(enml));
    assert(nrm > 0.99 && nrm < 1.01);
  }
#endif

  // Persistent storage.
  vtx_in_cell_lvl.reset(nlvl, t.nvtx);

  // Edge normal velocity -> vertex velocity.
  v2cslotsz = (3*c2cslotsz)/2;
  v2c.reset(v2cslotsz, t.nvtx);
  v2ccnt.reset(t.nvtx);
  const int e2v_nedge = e2v_function == 0 ? t.vdegree : square(t.vdegree);
  e2v_stencils.reset(e2v_nedge, t.nvtx);
  e2v_coefs.reset(3*e2v_nedge, t.nvtx);
  ompparfor for (int vi = 0; vi < t.nvtx; ++vi) {
    {
      int* vi_v2c = v2c(vi);
      std::set<int> unique;
      int cnt = 0;
      // The first three are the vertex's adjacent cells.
      bool halo_edge = false;
      for (int i = 0; i < t.vdegree; ++i) {
        const int ci = t.get_v2c(vi,i);
        if (ci < 0) {
          halo_edge = true;
          break;
        }
        unique.insert(ci);
        vi_v2c[cnt++] = ci;
      }
      if (halo_edge) {
        // No need to collect data for this vertex.
        v2ccnt[vi] = 0;
        continue; // vi loop
      }
      for (int i = 0; i < t.vdegree; ++i) {
        const int ci = t.get_v2c(vi,i);
        if (ci < 0) continue;
        const int nc2c = c2ccnt[ci];
        const int* ci_c2c = c2c(ci);
        for (int j = 0; j < nc2c; ++j) {
          const int cj = ci_c2c[j];
          if (cj < 0) continue; // halo edge
          if (unique.find(cj) != unique.end()) continue;
          unique.insert(cj);
          vi_v2c[cnt++] = cj;
        }
      }
      assert(cnt <= v2cslotsz);
      v2ccnt[vi] = cnt;
    }

    Vec xaxis, yaxis;
    Real xs[max_vdegree*max_vdegree], ys[max_vdegree*max_vdegree];
    int* stencil = e2v_stencils(vi);
    map_edges_to_vertex_plane(*mm, vi, xaxis, yaxis, xs, ys, stencil,
                              e2v_function == 0 ? 1 : 2);

    // Least-squares operator to reconstruct velocity at vertex.
    Real* B = e2v_coefs(vi);
    {
      static const int max_nrow = max_vdegree*max_vdegree + 3;
      static const int max_ncol = 9;
      const int nrow = e2v_nedge + (e2v_function == 0 ? 1 : 3);
      const int ncol = e2v_function == 0 ? 3 : 9;
      assert(nrow <= max_nrow);
      Vec vtx; fill(mm->v, vi, vtx, true);
      Real A[max_nrow*max_ncol] = {0};
      for (int i = 0; i < e2v_nedge; ++i) {
        assert(stencil[i] >= 0);
        const Real* enml = &enmls[3*stencil[i]];
        // The row encodes
        //     nml'(v0 + x vx + y vy) = enml,
        // Rearranging,
        //     [nml', x nml', y nml'] [v0; vx; vy] = enml,
        // and [v0; vx; vy] is the vector of unknowns. If e2v_function = 0, then
        // this is instead just nml'v0 = enml.
        //   It's not necessary to bring enml into the tangent plane, but it
        // makes sense to do so. Then we normalize it so the full edge normal
        // velocity is fit.
        Vec tenml;
        if (i >= t.vdegree) {
          geo::axpbyz(-geo::dot(enml, vtx), vtx, 1, enml, tenml);
          geo::normalize(tenml);
        } else {
          // The normals to the edges at this vtx are already tangent.
          geo::copy(tenml, enml);
        }
        for (int d = 0; d < 3; ++d) A[nrow*d + i] = tenml[d];
        if (e2v_function > 0) {
          for (int d = 0; d < 3; ++d) A[nrow*(d+3) + i] = xs[i]*tenml[d];
          for (int d = 0; d < 3; ++d) A[nrow*(d+6) + i] = ys[i]*tenml[d];
        }
      }
      // Add the tangency constraint as additional rows.
      for (int d = 0; d < 3; ++d) A[nrow*d + e2v_nedge] = vtx[d];
      if (e2v_function > 0) {
        // Each vector in v0 + x vx + y vy should be normal to vtx.
        for (int d = 0; d < 3; ++d) A[nrow*(d+3) + e2v_nedge+1] = vtx[d];
        for (int d = 0; d < 3; ++d) A[nrow*(d+6) + e2v_nedge+2] = vtx[d];
      }
      Real Btmp[max_nrow*max_ncol] = {0};
      if (e2v_function == 1) {
        // Solve as an unconstrained LS problem.
        Real tau[max_nrow*max_ncol], wrk[max_nrow], Q[max_nrow*max_ncol];
        lsinv(A, nrow, ncol, tau, wrk, Q, Btmp);
        // Since the last component is always 0 (vel vector normal to sphere
        // vtx), we can ignore the last col.
        //   If e2v_function > 0, we're interested in just the v0 part, so we
        // also can ignore rows 4 and up.
        for (int d = 0; d < 3; ++d)
          for (int i = 0; i < e2v_nedge; ++i)
            B[3*i+d] = Btmp[nrow*d+i];
      } else {
        // Solve as a constrained LS problem. The equality constraint duplicates
        // the final row of A. We keep this row because the upper square block
        // is rank deficient because the normals of the edges that meet at this
        // vertex are all in the tangent plane.
        assert(max_vdegree == 3);
        assert(t.vdegree == 3);
        for (int i = 0; i < nrow; ++i)
          Btmp[nrow*i+i] = 1;
        Real c[max_ncol];
        //todo To support e2v_function, need to permit an mx3 constraint matrix.
        for (int i = 0; i < ncol; ++i)
          c[i] = A[nrow*i + e2v_nedge];
        Real dy[max_ncol] = {0}, rwrk[2*max_nrow*max_ncol];
        solve_ls_1eq(nrow, ncol, ncol, A, c, Btmp, dy, rwrk);
        for (int i = 0; i < e2v_nedge; ++i)
          for (int d = 0; d < 3; ++d)
            B[3*i+d] = Btmp[nrow*i + d];
      }
#ifndef NDEBUG
      {
        const Real tol = ((e2v_function == 0 ? 10 : 100)*
                          std::numeric_limits<Real>::epsilon());
        for (int i = 0; i < 3; ++i) {
          Real a = 0;
          for (int d = 0; d < 3; ++d) a += B[3*i+d]*vtx[d];
          if (!(a <= tol)) prc(a);
          assert(a <= tol);
        }
      }
#endif
    }
  }

  // Interpolate vertex velocities to cell interior.
  cvsf.init(mm, he, 2, 8, false /* owned and shared */);
}
