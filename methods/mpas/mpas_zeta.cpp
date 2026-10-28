#include "mpas_zeta.hpp"
#include "mpas_geometry.hpp"

ZetaRep::ZetaRep (const MpasMesh::ConstPtr& mm, const siqp::sh::Mesh& m,
                  const Comm::Ptr& he) {
  const auto& t = mm->t;
  const auto& l = mm->l;
  const int nilvl = l.nmax+1;

  const auto& zint_levels = l.get_zint_array();
  if (zint_levels.m() == nilvl && zint_levels.n() == t.ncell)
    eul_cell = zint_levels;
  else {
    Ar2 z(nilvl, t.ncell);
    ompparfor for (int ci = 0; ci < t.ncell; ++ci) {
      const Real* zint = l.get_interface_z(ci);
      for (int li = 0; li < nilvl; ++li)
        z(li,ci) = zint[li];
    }
    eul_cell = z;
  }
  
  Ar3 zeta_eul_cell(1, nilvl, t.ncell);
  ompfor for (int ci = 0; ci < t.ncell; ++ci) {
    const Real* zint = l.get_interface_z(ci);
    for (int li = 0; li < nilvl; ++li)
      zeta_eul_cell(0,li,ci) = zint[li];
  }

  // We can move these back to being class data if we need these in Trajectory.
  CellFromCellsModalShapeFunction ccsf_eul;
  Ar3 shapes_eul;

  // Order 2 is consistent with the reconstruction in LagrangianLevels. In
  // particular, if dt = 0, then the linear LS fit in LL will recover the
  // cell-center Eulerian value. A limiter is permitted since that doesn't
  // change the mean.
  const bool use_fv_limiter = true;
  ccsf_eul.init(mm, m, he, 2, 2 /* unused */, use_fv_limiter, false);
  shapes_eul.reset(ccsf_eul.ho_ncoef, nilvl, t.ncell);
  eul_vtx_multivalued.reset(max_vtx, nilvl, t.ncell);
  ompfor for (int ci = 0; ci < t.ncell; ++ci) {
    if (at_halo_edge(*he, ci)) continue;
    for (int li = 0; li < nilvl; ++li) {
      Real* const shape = shapes_eul(li,ci);
      ccsf_eul.calc_shape(ci, zeta_eul_cell, 0, li, shape);
      const int ne = t.c2ne[ci];
      Real zmin, zmax;
      if (use_fv_limiter) {
        const Real zeta_lo = zeta_eul_cell(0,li,ci);
        zmin = zeta_lo;
        zmax = zmin;
        for (int j = 0; j < ne; ++j) {
          const int cj = t.get_c2c(ci,j);
          const Real zj = zeta_eul_cell(0,li,cj);
          zmin = std::min(zmin, zj);
          zmax = std::max(zmax, zj);
        }
        Real alpha = 1;
        const Real* const cc = ccsf_eul.ho_coefs(ci);
        const Real* const xaxis = cc + 3;
        const Real* const yaxis = xaxis + 3;
        for (int j = 0; j < ne; ++j) {
          const auto v = m.p(m.e(j,ci));
          Vec vmcc; geo::axpbyz(1, v, -1, cc, vmcc);
          const Real
            x = geo::dot(vmcc, xaxis),
            y = geo::dot(vmcc, yaxis),
            zeta_hi = zeta_lo + shape[0]*x + shape[1]*y;
          const bool above = zeta_hi > zmax, below = zeta_hi < zmin;
          if (below || above) {
            Real aj = (zmin == zmax ?
                       0 :
                       ((below ? zmin : zmax) - zeta_lo)/
                       (zeta_hi - zeta_lo));
            assert(aj <= 1);
            if (aj < 0) aj = 0;
            alpha = std::min(alpha, aj);
          }
        }
        assert(alpha >= 0 && alpha <= 1);
        if (alpha < 1)
          for (int j = 0; j < 2; ++j)
            shape[j] = alpha*shape[j];
      }
      ccsf_eul.eval_shape_at_vtxs(ci, zeta_eul_cell(0,li,ci), shape, 0, li,
                                  eul_vtx_multivalued(li,ci));
      if (use_fv_limiter) {
        for (int j = 0; j < ne; ++j)
          assert(eul_vtx_multivalued(j,li,ci) >= (1 - 1e-10)*zmin &&
                 eul_vtx_multivalued(j,li,ci) <= (1 + 1e-10)*zmax);
      }
    }
  }

#if 0
  // The univalued vertex value is the mean of the adjacent cell-center values.
  eul_vtx_univalued.reset(nilvl, t.nvtx);
  ompparfor for (int vi = 0; vi < t.nvtx; ++vi)
    fill_ent_zeta_eul(0, *mm, vi, eul_vtx_univalued(vi));
#else
  // The univalued vertex value is the mean of the multivalued vertex values.
  eul_vtx_univalued.reset(nilvl, t.nvtx);
  eul_vtx_univalued.set(0);
  ompparfor for (int vi = 0; vi < t.nvtx; ++vi) {
    int n = 0;
    for (int i = 0; i < t.vdegree; ++i) {
      const int ci = t.get_v2c(vi,i);
      if (ci < 0) continue;
      const int ne = t.c2ne[ci];
      bool found = false;
      for (int j = 0; j < ne; ++j) {
        if (m.e(j,ci) != vi) continue;
        found = true;
        ++n;
        for (int li = 0; li < nilvl; ++li)
          eul_vtx_univalued(li,vi) += eul_vtx_multivalued(j,li,ci);
        break;
      }
      assert(found);
    }
    assert(n > 0 && n <= max_vdegree);
    for (int li = 0; li < nilvl; ++li)
      eul_vtx_univalued(li,vi) /= n;
  }
#endif
}
