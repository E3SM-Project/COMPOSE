#include "mpas_mesh.hpp"
#include "mpas_geometry.hpp"
#include "mpas_zeta.hpp"
#include "mpas_comm.hpp"
#include "mpas_gallery.hpp"

#include <set>

MpasMesh
::MpasMesh (int ncell, int nedge, int nvtx, int max_c2e, int vdegree,
            const PyInt* c2ne, const PyInt* c2c, const PyInt* c2e, const PyInt* c2v,
            const PyInt* e2c, const PyInt* e2v,
            const PyInt* v2c, const PyInt* v2e,
            const Real* xc, const Real* yc, const Real* zc,
            const Real* xe, const Real* ye, const Real* ze,
            const Real* xv, const Real* yv, const Real* zv,
            const Real* dvEdge, const Real* dcEdge, const Real* angleEdge, const Real* areaCell,
            // For edge velocity reconstruction
            const PyInt* e2e, const PyInt* e2ne, const Real* e2w, int max_e2w,
            const int nlevel)
: t(ncell, nedge, nvtx, max_c2e, vdegree, c2ne, c2c, c2e, c2v, e2c, e2v, v2c, v2e,
    e2e, e2ne, e2w, max_e2w),
  c(ncell, xc, yc, zc),
  e(nedge, xe, ye, ze),
  v(nvtx, xv, yv, zv),
  g(dvEdge, dcEdge, angleEdge, areaCell),
  l(ncell, nedge, nlevel)
{
  // default these:
  ncown = ncell; neown = nedge; nvown = nvtx;
}

Levels::Levels (const int ncell, const int nedge, const int nlevel) {
  nmax = nlevel;
  dmin.resize(ncell + 2*nedge, 1);
  dmax.resize(ncell + 2*nedge, nlevel);
  cmin = dmin.data(); cmax = dmax.data();
  eminbot = cmin + ncell; emaxbot = cmax + ncell;
  emintop = eminbot + nedge; emaxtop = emaxbot + nedge;
}

void Levels::set_interface_z_all_cols (const Real* zint_) {
  Ar2 z(nmax+1, 1);
  for (int i = 0; i <= nmax; ++i) z(i,0) = zint_[i];
  zint = z;
}

void Levels
::set_interface_spatially_varying_zeta_eul (const CellCtrs& c, const Real zlo,
                                            const Real zhi) {
  Ar2 z(nmax+1, c.n);
  ompparfor for (int ci = 0; ci < c.n; ++ci) {
    const Real ctr = 0.75*(std::sin(M_PI*c.x[ci])*std::sin(M_PI*c.y[ci])*
                           std::sin(M_PI*c.z[ci]));
    for (int li = 0; li <= nmax; ++li) {
      const Real x = 2*(Real(li)/nmax) - 1;
      const Real f = get_level_shape_fn(x, ctr);
      z(li,ci) = zlo + (zhi - zlo)*f;
    }
  }
  zint = z;
}

const Real* Levels::get_interface_z (int c) const {
  if (zint.empty()) return nullptr;
  if (zint.n() == 1) return zint(0);
  assert(c >= 0 && c < zint.n());
  return zint(c);
}

void get_oriented_cell (const Topo& t, const Vtxs& v, const int ci, int* vtxs, int* edges) {
  const int ne = t.c2ne[ci];
  for (int i = 0; i < ne; ++i) {
    const int ei = t.get_c2e(ci,i), vi0 = t.get_e2v(ei,0), vi1 = t.get_e2v(ei,1);
    const int ej = t.get_c2e(ci, (i+1)%ne), vj0 = t.get_e2v(ej,0), vj1 = t.get_e2v(ej,1);
    assert(vi0 == vj0 || vi0 == vj1 || vi1 == vj0 || vi1 == vj1);
    vtxs[i] = vi1 == vj0 || vi1 == vj1 ? vi0 : vi1;
    if (edges) edges[i] = ei;
  }
}

void init_extras (MpasMesh& m) {
  const auto& t = m.t;
  auto& g = m.g;
  g.my_angleEdge.resize(t.nedge);
  g.my_sign.resize(t.ncell * t.max_vtx);
  g.my_areaCell.resize(t.ncell);
  g.my_cellCentroid.reset(3, t.ncell);

  Vec axis = {0,0,1};
  ompparfor for (int ei = 0; ei < t.nedge; ++ei) {
    Vec ec; fill(m.e, ei, ec, true);

    const int vi0 = t.get_e2v(ei,0), vi1 = t.get_e2v(ei,1);
    Vec v0; fill(m.v, vi0, v0, true);
    Vec v1; fill(m.v, vi1, v1, true);
    Vec vdif; geo::axpbyz(1, v1, -1, v0, vdif);
    Vec nml; geo::cross(vdif, ec, nml); geo::normalize(nml);

    Vec ec2; ec2[0] = ec[0]; ec2[1] = ec[1]; ec2[2] = 0;
    Vec east; geo::cross(axis, ec2, east); geo::normalize(east);
    Vec north; geo::cross(ec, east, north); geo::normalize(north);

    g.my_angleEdge[ei] = std::atan2(geo::dot(nml, north), geo::dot(nml, east));
  }

  ompparfor for (int ci = 0; ci < t.ncell; ++ci) {
    Vec cc; fill(m.c, ci, cc, true);
    const int ne = m.t.c2ne[ci];
    for (int i = 0; i < ne; ++i) {
      const int ei = m.t.get_c2e(ci,i);
      Vec ec; fill(m.e, ei, ec, true);

      const int vi0 = t.get_e2v(ei,0), vi1 = t.get_e2v(ei,1);
      Vec v0; fill(m.v, vi0, v0, true);
      Vec v1; fill(m.v, vi1, v1, true);
      Vec vdif; geo::axpbyz(1, v1, -1, v0, vdif);
      Vec enml; geo::cross(vdif, ec, enml); geo::normalize(enml);

      Vec ecmcc; geo::axpbyz(1, ec, -1, cc, ecmcc); geo::normalize(ecmcc);

      const short sign = geo::dot(enml, ecmcc) >= 0 ? 1 : -1;
      g.my_sign[ci*t.max_vtx + i] = sign;
    }
  }
}

siqp::sh::Mesh::Ptr create_siqp_mesh (const MpasMesh& mm) {
  const auto& t = mm.t;
  const auto& v = mm.v;
  Ar2 p(3, t.nvtx);
  Array2D<int> e(max_vtx, t.ncell);
  omppar {
    ompfor for (int i = 0; i < t.nvtx; ++i) {
      p(0,i) = v.x[i]; p(1,i) = v.y[i]; p(2,i) = v.z[i];
      geo::normalize(p(i));
    }
    ompfor for (int ci = 0; ci < t.ncell; ++ci) {
      const int ne = t.c2ne[ci];
      assert(ne <= max_vtx);
      for (int i = ne; i < max_vtx; ++i) e(i,ci) = -1;
      get_oriented_cell(t, v, ci, e(ci));
    }
  }
  auto m = std::make_shared<siqp::sh::Mesh>();
  m->p = p;
  m->e = e;
  siqp::fill_normals(*m);
  return m;
}

void init_more_extras (MpasMesh& mm, siqp::sh::Mesh& m) {
  const auto& t = mm.t;
  mm.g.my_areaCell.resize(t.ncell);
  mm.g.my_cellCentroid.reset(3, t.ncell);
  const TriangleQuadrature tq;
  ompparfor for (int ci = 0; ci < t.ncell; ++ci) {
    const int ne = t.c2ne[ci];
    assert(ne <= 20);
    Real vs[3*20];
    for (int i = 0; i < ne; ++i) geo::copy(vs+3*i, m.p(m.e(i,ci)));
    mm.g.my_areaCell[ci] = geo::calc_area(tq, vs, ne);
    calc_cell_centroid(mm, ci, slice(mm.g.my_cellCentroid, ci));
  }
}

void init_zeta_rep (const MpasMesh::Ptr& mm, const siqp::sh::Mesh& m,
                    const Comm::Ptr& he) {
  mm->zeta_rep = std::make_shared<ZetaRep>(mm, m, he);
}

void calc_tangent_plane_axes (Vec ctr, Vec xaxis, Vec yaxis) {
  Vec up = {0,0,1}, right = {1,0,0};
  if (std::abs(ctr[2]) > 0.707) {
    geo::cross(right, ctr, xaxis); geo::normalize(xaxis);
    geo::cross(ctr, xaxis, yaxis); geo::normalize(yaxis);
  } else {
    geo::cross(up, ctr, xaxis); geo::normalize(xaxis);
    geo::cross(ctr, xaxis, yaxis); geo::normalize(yaxis);
  }
}

void calc_xy (const Vec ctr, const Vec p, const Vec xaxis, const Vec yaxis,
              const bool scale_by_arc_length, Real& x, Real& y) {
  Vec pmctr; geo::axpbyz(1, p, -1, ctr, pmctr);
  x = geo::dot(pmctr, xaxis);
  y = geo::dot(pmctr, yaxis);
  if (scale_by_arc_length) {
    // Scale by arc length.
    const Real scale = (geo::calc_arc_length(p, ctr) /
                        std::sqrt(square(x) + square(y)));
    x *= scale;
    y *= scale;
  }
}

int map_cells_to_cell_plane (const MpasMesh& m, const int ci,
                             Vec xaxis, Vec yaxis, Real* xs, Real* ys,
                             int* stencil, const bool scale_by_arc_length,
                             const bool two_halo) {
  const auto& t = m.t;
  const auto& c = m.c;

  // Cell-center local x-y coordinate system.
  Vec cc; fill(c, ci, cc, true);
  calc_tangent_plane_axes(cc, xaxis, yaxis);

  // Project neighbor cell centers onto tangent plane.
  xs[0] = 0; ys[0] = 0; // central point is always at the origin
  const int ne = t.c2ne[ci];
  for (int i = 0; i < ne; ++i) {
    const int ei = t.get_c2e(ci,i);
    const int cother = get_other_cell(t, ci, ei);
    assert(cother >= 0);
    if (stencil) stencil[i] = cother;
    Vec co; fill(c, cother, co, true);
    calc_xy(cc, co, xaxis, yaxis, scale_by_arc_length, xs[i+1], ys[i+1]);
  }

  if ( ! two_halo) return ne;

  assert(stencil);
  std::set<int> unique;
  unique.insert(ci);
  for (int i = 0; i < ne; ++i) unique.insert(stencil[i]);
  int cnt = ne;
  for (int i = 0; i < ne; ++i) {
    const int ei = t.get_c2e(ci,i);
    const int cother = get_other_cell(t, ci, ei);
    assert(cother >= 0);
    const int one = t.c2ne[cother];
    for (int j = 0; j < one; ++j) {
      const int ej = t.get_c2e(cother,j);
      const int cj = get_other_cell(t, cother, ej);
      assert(cj >= 0);
      if (unique.find(cj) != unique.end()) continue;
      unique.insert(cj);
      assert(cnt < max_2halo_slots-1);
      stencil[cnt] = cj;
      Vec ccj; fill(c, cj, ccj, true);
      calc_xy(cc, ccj, xaxis, yaxis, scale_by_arc_length, xs[cnt+1], ys[cnt+1]);
      ++cnt;
    }
  }
  return cnt;
}

void map_vtxs_to_cell_plane (
  const MpasMesh& m, const int ci, Vec xaxis, Vec yaxis, Real* xs, Real* ys,
  int* stencil, const bool scale_by_arc_length, const Real* cc_replacement,
  const Array3D<const Real>* v_replacement, const int lvl)
{
  assert((cc_replacement && v_replacement) || ! (cc_replacement && v_replacement));

  const auto& t = m.t;
  const auto& c = m.c;
  const auto& v = m.v;

  // Cell-center local x-y coordinate system.
  Vec cc;
  if (cc_replacement) geo::copy(cc, cc_replacement);
  else fill(c, ci, cc, true);
  calc_tangent_plane_axes(cc, xaxis, yaxis);

  // Project vertices onto tangent plane.
  const int ne = t.c2ne[ci];
  for (int i = 0; i < ne; ++i) {
    const int vi = t.get_c2v(ci,i);
    Vec vtx;
    if (v_replacement) geo::copy(vtx, (*v_replacement)(lvl,vi));
    else fill(v, vi, vtx, true);
    if (stencil) stencil[i] = vi;
    calc_xy(cc, vtx, xaxis, yaxis, scale_by_arc_length, xs[i], ys[i]);
  }  
}

void map_edges_to_vertex_plane (const MpasMesh& m, const int vi,
                                Vec xaxis, Vec yaxis, Real* xs, Real* ys,
                                int* stencil, const int nhop,
                                const bool scale_by_arc_length) {
  const auto& t = m.t;
  const auto& v = m.v;
  const auto& e = m.e;

  // Vertex-local x-y coordinate system.
  Vec vtx; fill(v, vi, vtx, true);
  calc_tangent_plane_axes(vtx, xaxis, yaxis);

  // Project neighbor edge centers onto tangent plane.
  for (int i = 0; i < t.vdegree; ++i) {
    const int ei = t.get_v2e(vi,i);
    assert(ei >= 0);
    if (stencil) stencil[i] = ei;
    Vec ec; fill(e, ei, ec, true);
    calc_xy(vtx, ec, xaxis, yaxis, scale_by_arc_length, xs[i], ys[i]);
  }
  if (nhop > 1) {
    throw_if(nhop > 2, "only nhop <= 2 supported");
    throw_if( ! stencil, "if nhop > 1, then stencil must be non-null");
    int idx = t.vdegree;
    for (int j = 0; j < t.vdegree; ++j) {
      const int ej = stencil[j];
      const int vj = t.get_e2v(ej,0) == vi ? t.get_e2v(ej,1) : t.get_e2v(ej,0);
      assert(vj >= 0);
      for (int k = 0; k < t.vdegree; ++k) {
        const int ek = t.get_v2e(vj,k);
        if (ek == ej) continue;
        stencil[idx] = ek;
        Vec ec; fill(e, ek, ec, true);
        calc_xy(vtx, ec, xaxis, yaxis, scale_by_arc_length, xs[idx], ys[idx]);
        ++idx;
      }
    }
    assert(idx == square(t.vdegree));
  }
}

void calc_cell_centroid (const MpasMesh& m, const int ci, Vec cc) {
  const int ne = m.t.c2ne[ci];

  int elem[max_vtx];
  get_oriented_cell(m.t, m.v, ci, elem);

  Real vs[3*max_vtx];
  for (int i = 0; i < ne; ++i) fill(m.v, elem[i], vs + 3*i, true);

  TriangleQuadrature tq;
  const Real* qx, * qw;
  int nq;
  tq.get_coef(12, qx, qw, nq);

  Vec mcc;
  fill(m.c, ci, mcc, true);

  geo::set(cc, 0);
  for (int i = 0; i < ne; ++i) {
    for (Int k = 0; k < nq; ++k) {
      Vec u;
      const Real jac = geo::calc_tri_jacobian(
        mcc, slice(vs,i), slice(vs, (i+1) % ne),
        slice(qx, k), u);
      for (int d = 0; d < 3; ++d) cc[d] += qw[k]*jac*u[d];
    }
  }
  geo::normalize(cc);
}

int augment_stencil (const MpasMesh& m, const int ci, int ncurrent, int* stencil,
                     const int nhalo, const bool permit_neg_cid) {
  assert(stencil);
  std::set<int> unique;
  unique.insert(ci);
  for (int i = 0; i < ncurrent; ++i) unique.insert(stencil[i]);
  int start = ncurrent == 0 ? -1 : 0, cnt = ncurrent;
  for (int halo = 0; halo < nhalo; ++halo) {
    for (int i = start; i < ncurrent; ++i) {
      const int cother = start == -1 ? ci : stencil[i];
      const int one = m.t.c2ne[cother];
      for (int j = 0; j < one; ++j) {
        const int ej = m.t.get_c2e(cother,j);
        const int cj = get_other_cell(m.t, cother, ej);
        if (permit_neg_cid && cj < 0) continue;
        if (unique.find(cj) != unique.end()) continue;
        unique.insert(cj);
        stencil[cnt++] = cj;
      }
    }
    start = ncurrent;
    ncurrent = cnt;
  }
  return cnt;
}
