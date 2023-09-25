#include "slmm_mesh.hpp"
#include "slmm_gll.hpp"
#include "slmm_util.hpp"
#include "slmm_basis.hpp"
#include "slmm_islet.hpp"
using geometry = siqk::SphereGeometry;

#include <map>
#include <set>
#include <stdexcept>

namespace slmm {
namespace mesh {

static void make_equiangular_nodes (const Int ne, std::vector<Real>& x) {
  static const Real oosqrt3 = 1.0 / std::sqrt(3.0);
  const Real dtheta = 0.5*M_PI / ne;
  x.resize(ne+1);
  if (ne % 2 == 1) {
    const Int n = (ne + 1) / 2;
    for (Int i = 0; i < n; ++i)
      x[n + i] = oosqrt3*std::tan((i + 0.5)*dtheta);
    for (Int i = 0; i < n; ++i)
      x[n - 1 - i] = -x[n + i];
  } else {
    const Int n = ne / 2;
    x[n] = 0;
    for (Int i = 1; i <= n; ++i)
      x[n + i] = oosqrt3*std::tan(i*dtheta);
    for (Int i = 1; i <= n; ++i)
      x[n - i] = -x[n + i];
  }
}

static void make_planar_mesh (AVec3s& p, AIdxs& e, const Int n) {
  std::vector<Real> x;
  make_equiangular_nodes(n, x);
  resize(e, n*n, 4);
  resize(p, (n+1)*(n+1));
  for (Int iy = 0; iy < n+1; ++iy)
    for (Int ix = 0; ix < n+1; ++ix) {
      const auto idx = (n+1)*iy + ix;
      p(idx,0) = x[ix];
      p(idx,1) = x[iy];
      p(idx,2) = 0;
    }
  for (Int iy = 0; iy < n; ++iy)
    for (Int ix = 0; ix < n; ++ix) {
      const auto idx = n*iy + ix;
      e(idx,0) = (n+1)*iy + ix;
      e(idx,1) = (n+1)*iy + ix+1;
      e(idx,2) = (n+1)*(iy+1) + ix+1;
      e(idx,3) = (n+1)*(iy+1) + ix;
    }
}

template <typename V>
static void rotate (const Real R[9], V p) {
  const Real x = p[0], y = p[1], z = p[2];
  p[0] = R[0]*x + R[1]*y + R[2]*z;
  p[1] = R[3]*x + R[4]*y + R[5]*z;
  p[2] = R[6]*x + R[7]*y + R[8]*z;
}

template <typename V>
static void translate (const Real xlate[3], V p) {
  for (Int i = 0; i < 3; ++i) p[i] += xlate[i];
}

static void transform_planar_mesh (const Real R[9], const Real xlate[3],
                                   AVec3s& p) {
  for (Int i = 0; i < nslices(p); ++i) {
    rotate(R, slice(p, i));
    translate(xlate, slice(p, i));
  }
}

// Remove vertices marked unused and adjust numbering.
static void remove_unused_vertices (AVec3s& p, AIdxs& e,
                                    const Real unused) {
  // adjust[i] is the number to subtract from i. Hence if e(ei,0) was originally
  // i, it is adjusted to i - adjust[i].
  std::vector<Int> adjust(nslices(p), 0);
  Int rmcnt = 0;
  for (Int i = 0; i < nslices(p); ++i) {
    if (p(i,0) != unused) continue;
    adjust[i] = 1;
    ++rmcnt;
  }
  // Cumsum.
  for (Int i = 1; i < nslices(p); ++i)
    adjust[i] += adjust[i-1];
  // Adjust e.
  for (Int ei = 0; ei < nslices(e); ++ei)
    for (Int k = 0; k < szslice(e); ++k)
      e(ei,k) -= adjust[e(ei,k)];
  // Remove unused from p.
  AVec3s pc(nslices(p));
  copy(pc, p);
  resize(p, nslices(p) - rmcnt);
  for (Int i = 0, j = 0; i < nslices(pc); ++i) {
    if (pc(i,0) == unused) continue;
    for (Int k = 0; k < szslice(pc); ++k) p(j,k) = pc(i,k);
    ++j;
  }
}

void make_cubedsphere_mesh (AVec3s& p, AIdxs& e, const Int n) {
  // Transformation of the reference mesh make_planar_mesh to make each of the
  // six faces.
  const Real d = 1.0 / std::sqrt(3.0);
  static Real R[6][9] = {{ 1, 0, 0, 0, 0, 0, 0, 1, 0},  // face 0, -y
                         { 0, 0, 0, 1, 0, 0, 0, 1, 0},  //      1, +x
                         {-1, 0, 0, 0, 0, 0, 0, 1, 0},  //      2, +y
                         { 0, 0, 0,-1, 0, 0, 0, 1, 0},  //      3, -x
                         { 1, 0, 0, 0, 1, 0, 0, 0, 0},  //      4, +z
                         {-1, 0, 0, 0, 1, 0, 0, 0, 0}}; //      5, -z
  static Real xlate[6][3] = {{ 0,-d, 0}, { d, 0, 0}, { 0, d, 0},
                             {-d, 0, 0}, { 0, 0, d}, { 0, 0,-d}};
  // Construct 6 uncoupled faces.
  AVec3s ps[6];
  AVec3s& p_ref = ps[0];
  AIdxs es[6];
  AIdxs& e_ref = es[0];
  make_planar_mesh(p_ref, e_ref, n);
  resize(e, 6*nslices(e_ref), 4);
  resize(p, 6*nslices(p_ref));
  for (Int i = 1; i < 6; ++i) {
    resize(es[i], nslices(e_ref), 4);
    copy(es[i], e_ref);
    resize(ps[i], nslices(p_ref));
    copy(ps[i], p_ref);
    transform_planar_mesh(R[i], xlate[i], ps[i]);
  }
  transform_planar_mesh(R[0], xlate[0], ps[0]);
  // Pack (p,e), accounting for equivalent vertices. For the moment, keep the p
  // slot for an equivalent vertex to make node numbering simpler, but make the
  // value bogus so we know if there's a problem in the numbering.
  const Real unused = -2;
  copy(p, unused);
  Int p_base = 0, e_base = 0;
  { // -y face
    const auto& fp = ps[0];
    auto& fe = es[0];
    for (Int j = 0; j < nslices(fp); ++j)
      for (Int k = 0; k < 3; ++k) p(j,k) = fp(j,k);
    for (Int j = 0; j < nslices(fe); ++j)
      for (Int k = 0; k < 4; ++k) e(j,k) = fe(j,k);
    p_base += nslices(p_ref);
    e_base += nslices(e_ref);
  }
  for (Int fi = 1; fi <= 2; ++fi) { // +x, +y faces
    const auto& fp = ps[fi];
    auto& fe = es[fi];
    for (Int j = 0; j < nslices(fp); ++j) {
      if (j % (n+1) == 0) continue; // equiv vertex
      for (Int k = 0; k < 3; ++k) p(p_base+j,k) = fp(j,k);
    }
    for (Int j = 0; j < nslices(fe); ++j) {
      for (Int k = 0; k < 4; ++k) fe(j,k) += p_base;
      // Left 2 vertices of left elem on face fi equiv to right 2 vertices of
      // right elem on face fi-1. Write to the face, then copy to e, so that
      // other faces can use these updated data.
      if (j % n == 0) {
        fe(j,0) = es[fi-1](j+n-1,1);
        fe(j,3) = es[fi-1](j+n-1,2);
      }
      for (Int k = 0; k < 4; ++k) e(e_base+j,k) = fe(j,k);
    }
    p_base += nslices(p_ref);
    e_base += nslices(e_ref);
  }
  { // -x face
    const auto& fp = ps[3];
    auto& fe = es[3];
    for (Int j = 0; j < nslices(fp); ++j) {
      if (j % (n+1) == 0 || (j+1) % (n+1) == 0) continue;
      for (Int k = 0; k < 3; ++k) p(p_base+j,k) = fp(j,k);
    }
    for (Int j = 0; j < nslices(fe); ++j) {
      for (Int k = 0; k < 4; ++k) fe(j,k) += p_base;
      if (j % n == 0) {
        fe(j,0) = es[2](j+n-1,1);
        fe(j,3) = es[2](j+n-1,2);
      } else if ((j+1) % n == 0) {
        fe(j,1) = es[0]((j+1)-n,0);
        fe(j,2) = es[0]((j+1)-n,3);
      }
      for (Int k = 0; k < 4; ++k) e(e_base+j,k) = fe(j,k);
    }
    p_base += nslices(p_ref);
    e_base += nslices(e_ref);
  }
  { // +z face
    const auto& fp = ps[4];
    auto& fe = es[4];
    for (Int j = n+1; j < nslices(fp) - (n+1); ++j) {
      if (j % (n+1) == 0 || (j+1) % (n+1) == 0) continue;
      for (Int k = 0; k < 3; ++k) p(p_base+j,k) = fp(j,k);
    }
    for (Int j = 0; j < nslices(fe); ++j)
      for (Int k = 0; k < 4; ++k) fe(j,k) += p_base;
    for (Int j = 0; j < n; ++j) { // -y
      fe(j,0) = es[0](n*(n-1)+j,3);
      fe(j,1) = es[0](n*(n-1)+j,2);
    }
    for (Int j = 0; j < n; ++j) { // +y
      fe(n*(n-1)+j,2) = es[2](n*n-1-j,3);
      fe(n*(n-1)+j,3) = es[2](n*n-1-j,2);
    }
    for (Int j = 0, i3 = 0; j < nslices(fe); j += n, ++i3) { // -x
      fe(j,0) = es[3](n*n-1-i3,2);
      fe(j,3) = es[3](n*n-1-i3,3);
    }
    for (Int j = n-1, i1 = 0; j < nslices(fe); j += n, ++i1) { // +x
      fe(j,1) = es[1](n*(n-1)+i1,3);
      fe(j,2) = es[1](n*(n-1)+i1,2);
    }
    for (Int j = 0; j < nslices(fe); ++j)
      for (Int k = 0; k < 4; ++k) e(e_base+j,k) = fe(j,k);
    p_base += nslices(p_ref);
    e_base += nslices(e_ref);
  }
  { // -z face
    const auto& fp = ps[5];
    auto& fe = es[5];
    for (Int j = n+1; j < nslices(fp) - (n+1); ++j) {
      if (j % (n+1) == 0 || (j+1) % (n+1) == 0) continue;
      for (Int k = 0; k < 3; ++k) p(p_base+j,k) = fp(j,k);
    }
    for (Int j = 0; j < nslices(fe); ++j)
      for (Int k = 0; k < 4; ++k) fe(j,k) += p_base;
    for (Int j = 0; j < n; ++j) { // -y
      fe(j,0) = es[0](n-1-j,1);
      fe(j,1) = es[0](n-1-j,0);
    }
    for (Int j = 0; j < n; ++j) { // +y
      fe(n*(n-1)+j,2) = es[2](j,1);
      fe(n*(n-1)+j,3) = es[2](j,0);
    }
    for (Int j = 0, i3 = 0; j < nslices(fe); j += n, ++i3) { // -x
      fe(j,0) = es[1](i3,0);
      fe(j,3) = es[1](i3,1);
    }
    for (Int j = n-1, i1 = 0; j < nslices(fe); j += n, ++i1) { // +x
      fe(j,1) = es[3](n-1-i1,1);
      fe(j,2) = es[3](n-1-i1,0);
    }
    for (Int j = 0; j < nslices(fe); ++j)
      for (Int k = 0; k < 4; ++k) e(e_base+j,k) = fe(j,k);
  }
  // Now go back and remove the unused vertices and adjust the numbering.
  remove_unused_vertices(p, e, unused);
  // Project to the unit sphere.
  for (Int i = 0; i < nslices(p); ++i)
    geometry::normalize(slice(p, i));
}

void make_subcell_from_geo (
  const AVec3s& geo_p, const AIdxs& geo_c2n, const Int np,
  const Real* const ref_x, AVec3s& sc_p, AIdxs& sc_c2n)
{
  AIdxs geo_c2e, geo_e2n;
  impl::make_c2e_from_c2n(np, geo_c2n, geo_c2e, geo_e2n);
  resize(sc_p,
         nslices(geo_p) +                      // corner nodes
         (np-2)*nslices(geo_e2n) +             // np-2 per edge
         siqk::square(np-2)*nslices(geo_c2n)); // nodes inside a geo cell
  resize(sc_c2n, nslices(geo_c2n), siqk::square(np));
  Int pi = 0;

  // Geo cell vertices.
  for ( ; pi < nslices(geo_p); ++pi)
    for (Int k = 0; k < 3; ++k)
      sc_p(pi,k) = geo_p(pi,k);
  ko::View<int**, ko::HostSpace> nodes("nodes", np, np);

  // Add new edge nodes.
  for (Int gci = 0; gci < nslices(geo_c2n); ++gci) {
    const auto geo_nodes = slice(geo_c2n, gci);
    for (Int ei = 0; ei < 4; ++ei) {
      // If my edge is i -> j and j > i, then I'm responsible for adding these.
      if (geo_nodes[ei] > geo_nodes[(ei+1) % 4]) continue;
      // edge[0] -> edge[np-1] is the geo edge.
      auto edge = slice(geo_e2n, geo_c2e(gci,ei));
      assert(edge[0] == geo_nodes[ei]);
      assert(edge[np-1] == geo_nodes[(ei+1) % 4]);
      // Add the new nodes.
      const auto p0 = slice(sc_p, edge[0]);
      const auto p1 = slice(sc_p, edge[np-1]);
      for (Int i = 1; i < np-1; ++i) {
        auto p = slice(sc_p, pi);
        const Real alpha = 0.5*(ref_x[i] + 1);
        for (Int k = 0; k < 3; ++k)
          p[k] = (1 - alpha)*p0[k] + alpha*p1[k];
        edge[i] = pi;
        ++pi; 
      }
    }
  }
  for (Int gci = 0; gci < nslices(geo_c2n); ++gci) {
    const auto geo_nodes = slice(geo_c2n, gci);
    // Record the newly created edge nodes.
    for (Int ei = 0; ei < 4; ++ei) {
      const auto edge = slice(geo_e2n, geo_c2e(gci,ei));
      if (geo_nodes[ei] < geo_nodes[(ei+1) % 4]) {
        assert(edge[0] == geo_nodes[ei]);
        assert(edge[np-1] == geo_nodes[(ei+1) % 4]);
        switch (ei) {
        case 0: for (Int i = 0; i < np; ++i) nodes(i,0)    = edge[i];      break;
        case 1: for (Int i = 0; i < np; ++i) nodes(np-1,i) = edge[i];      break;
        case 2: for (Int i = 0; i < np; ++i) nodes(i,np-1) = edge[np-1-i]; break;
        case 3: for (Int i = 0; i < np; ++i) nodes(0,i)    = edge[np-1-i]; break;
        default: assert(0);
        }
      } else {
        assert(edge[np-1] == geo_nodes[ei]);
        assert(edge[0] == geo_nodes[(ei+1) % 4]);
        switch (ei) {
        case 0: for (Int i = 0; i < np; ++i) nodes(i,0)    = edge[np-1-i]; break;
        case 1: for (Int i = 0; i < np; ++i) nodes(np-1,i) = edge[np-1-i]; break;
        case 2: for (Int i = 0; i < np; ++i) nodes(i,np-1) = edge[i];      break;
        case 3: for (Int i = 0; i < np; ++i) nodes(0,i)    = edge[i];      break;
        default: assert(0);
        }
      }
    }
    // Add new internal nodes.
    for (Int j = 1; j < np-1; ++j) {
      const auto p0 = slice(sc_p, nodes(0,j));
      const auto p1 = slice(sc_p, nodes(np-1,j));
      for (Int i = 1; i < np-1; ++i) {
        assert(pi < nslices(sc_p));
        auto p = slice(sc_p, pi);
        const Real alpha = 0.5*(ref_x[i] + 1);
        for (Int k = 0; k < 3; ++k)
          p[k] = (1 - alpha)*p0[k] + alpha*p1[k];
        nodes(i,j) = pi;
        ++pi;
      }
    }
    // Fill CGLL cell with nodes.
    {
      auto cell = slice(sc_c2n, gci);
      for (Int j = 0, k = 0; j < np; ++j)
        for (Int i = 0; i < np; ++i, ++k)
          cell[k] = nodes(i,j);
    }
  }
  // Project to the unit sphere. (The first nodes are already projected.)
  for (Int i = nslices(geo_p); i < nslices(sc_p); ++i)
    geometry::normalize(slice(sc_p, i));
}  

static void make_subcell_from_geo (
  const AVec3s& geo_p, const AIdxs& geo_c2n, const Int np,
  const Basis& basis, AVec3s& sc_p, AIdxs& sc_c2n)
{
  // Reference cell nodes.
  const Real* ref_x = nullptr;
  ko::View<Real*, ko::HostSpace> buf("buf", np);
  basis.get_x(np, ref_x);
  make_subcell_from_geo(geo_p, geo_c2n, np, ref_x, sc_p, sc_c2n);
}

void make_cgll_from_geo (
  const AVec3s& geo_p, const AIdxs& geo_c2n, const Int np,
  const Basis& basis, AVec3s& cgll_p, AIdxs& cgll_c2n)
{
  make_subcell_from_geo(geo_p, geo_c2n, np, basis, cgll_p, cgll_c2n);
}

void make_io_cgll_from_internal_cgll (
  const AVec3s& cgll_p, const AIdxs& cgll_c2n,
  AIdxs& cgll_io_c2n)
{
  const Int np2 = szslice(cgll_c2n), np = std::sqrt(np2),
    nsc = siqk::square(np-1);
  resize(cgll_io_c2n, nslices(cgll_c2n)*nsc, 4);
  for (Int ci = 0; ci < nslices(cgll_c2n); ++ci) {
    const auto cell = slice(cgll_c2n, ci);
    for (Int scj = 0; scj < np-1; ++scj)
      for (Int sci = 0; sci < np-1; ++sci) {
        auto subcell = slice(cgll_io_c2n, nsc*ci + (np-1)*scj + sci);
        subcell[0] = cell[np* scj    + sci  ];
        subcell[1] = cell[np* scj    + sci+1];
        subcell[2] = cell[np*(scj+1) + sci+1];
        subcell[3] = cell[np*(scj+1) + sci  ];
      }
  }
}

void make_dgll_from_cgll (
  const AVec3s& cgll_p, const AIdxs& cgll_c2n,
  AIdxArray& dglln2cglln, AIdxs& dgll_c2n)
{
  const Int np2 = szslice(cgll_c2n);
  resize(dglln2cglln, np2*nslices(cgll_c2n));
  resize(dgll_c2n, nslices(cgll_c2n), szslice(cgll_c2n));
  AIdxArray cgll_c2n_used(nslices(cgll_p));
  copy(cgll_c2n_used, 0);
  Int pi = nslices(cgll_p);
  for (Int ci = 0; ci < nslices(cgll_c2n); ++ci) {
    const auto cgll_cell = slice(cgll_c2n, ci);
    auto dgll_cell = slice(dgll_c2n, ci);
    for (Int ni = 0; ni < np2; ++ni) {
      const Int cgll_node_nmbr = cgll_cell[ni];
      dglln2cglln[ci*np2 + ni] = cgll_node_nmbr;
      if (cgll_c2n_used[cgll_node_nmbr])
        dgll_cell[ni] = pi++;
      else {
        dgll_cell[ni] = cgll_node_nmbr;
        cgll_c2n_used[cgll_node_nmbr] = 1;
      }
    }
  }
#ifndef NDEBUG
  assert(pi == nslices(cgll_c2n) * np2);
  for (Int i = 0; i < nslices(cgll_c2n_used); ++i)
    assert(cgll_c2n_used[i]);
#endif
}

void make_cubedsphere_subcell_mesh (
  const Int ne, const Int np, const Basis& basis, AVec3s& geo_p,
  AIdxs& geo_c2n)
{
  AIdxs internal_c2n;
  {
    AVec3s cs_p;
    AIdxs cs_c2n;
    make_cubedsphere_mesh(cs_p, cs_c2n, ne);
    make_subcell_from_geo(cs_p, cs_c2n, np, basis, geo_p, internal_c2n);
  }
  make_io_cgll_from_internal_cgll(geo_p, internal_c2n, geo_c2n);
}

void get_adjacent_cells (
  const AIdxs& geo_c2n, AIdxArray& geo_c2cnbrs_ptr,
  AIdxArray& geo_c2cnbrs)
{
  const Int ncell = nslices(geo_c2n);
  std::map<Int,std::vector<Int> > n2cs;
  for (Int i = 0; i < ncell; ++i)
    for (Int k = 0; k < 4; ++k)
      n2cs[geo_c2n(i,k)].push_back(i);
  std::vector<std::set<Int> > c2cs(ncell);
  for (Int i = 0; i < ncell; ++i)
    for (Int k = 0; k < 4; ++k) {
      const auto& cells = n2cs[geo_c2n(i,k)];
      for (const auto cell : cells)
        c2cs[cell].insert(cells.begin(), cells.end());
    }
  {
    Int nn = 0;
    for (Int i = 0; i < ncell; ++i) {
      c2cs[i].erase(i);
      nn += c2cs[i].size();
    }
    resize(geo_c2cnbrs, nn);
  }
  resize(geo_c2cnbrs_ptr, ncell+1);
  geo_c2cnbrs_ptr[0] = 0;
  for (Int i = 0; i < ncell; ++i) {
    const Int nn = c2cs[i].size();
    geo_c2cnbrs_ptr[i+1] = geo_c2cnbrs_ptr[i] + nn;
    Int j = geo_c2cnbrs_ptr[i];
    for (const auto cell : c2cs[i])
      geo_c2cnbrs[j++] = cell;
  }
}

namespace impl {
void calc_elem_ctr (const AVec3s& p, const AIdxs& e,
                    const Int ei, Real ctr[3]) {
  for (Int j = 0; j < 3; ++j) ctr[j] = 0;
  Int n = 0;
  for (Int i = 0; i < szslice(e); ++i) {
    if (e(ei,i) < 0) break;
    for (Int j = 0; j < 3; ++j) ctr[j] += p(e(ei,i),j);
    ++n;
  }
  for (Int j = 0; j < 3; ++j) ctr[j] /= n;
}

struct Edge {
  const Int lo, hi;
  Edge (const Int& n0, const Int& n1)
    : lo(n0 < n1 ? n0 : n1),
      hi(n0 < n1 ? n1 : n0)
  {}
  bool operator< (const Edge& e) const {
    if (lo < e.lo) return true;
    if (lo == e.lo) return hi < e.hi;
    return false;
  }
};

void make_c2e_from_c2n (const Int np, const AIdxs& c2n,
                        AIdxs& c2e, AIdxs& e2n) {
  const Int nnode = szslice(c2n);
  // Number the edges.
  std::map<Edge, Int> edge2nmbr;
  Int nmbr = 0;
  for (Int ci = 0; ci < nslices(c2n); ++ci) {
    const auto cell = slice(c2n, ci);
    for (Int ni = 0; ni < nnode; ++ni) {
      Edge e(cell[ni], cell[(ni+1) % nnode]);
      const auto it = edge2nmbr.find(e);
      if (it == edge2nmbr.end())
        edge2nmbr[e] = nmbr++;
    }
  }
  // Fill the adjacency arrays.
  resize(c2e, nslices(c2n), szslice(c2n));
  resize(e2n, nmbr, np);
  for (Int ci = 0; ci < nslices(c2n); ++ci) {
    const auto cell = slice(c2n, ci);
    for (Int ni = 0; ni < nnode; ++ni) {
      Edge e(cell[ni], cell[(ni+1) % nnode]);
      const auto it = edge2nmbr.find(e);
      assert(it != edge2nmbr.end());
      const Int nmbr = it->second;
      c2e(ci, ni) = nmbr;
      e2n(nmbr, 0)    = it->first.lo;
      e2n(nmbr, np-1) = it->first.hi;
    }
  }
}

Int check_elem_normal_against_sphere (const AVec3s& p,
                                      const AIdxs& e) {
  Int nerr = 0;
  for (Int ei = 0; ei < nslices(e); ++ei) { // for each element
    Real sphere[3]; // ray through elem ctr
    calc_elem_ctr(p, e, ei, sphere);
    for (Int ti = 0; ti < szslice(e) - 2; ++ti) { // for each tri
      if (e(ei,ti+2) < 0) break;
      Real tri_normal[3]; {
        Real v[2][3];
        for (Int j = 0; j < 2; ++j) {
          geometry::copy(v[j], slice(p, e(ei,ti+j+1)));
          geometry::axpy(-1, slice(p, e(ei,0)), v[j]);
        }
        geometry::cross(v[0], v[1], tri_normal);
      }
      if (geometry::dot(tri_normal, sphere) <= 0)
        ++nerr;
    }
  }
  return nerr;
}
} // namespace impl

namespace impl {
// Split n things into two parts.
void split (const Int n, Int d[2]) {
  d[0] = n >> 1;
  d[1] = (d[0] << 1 == n) ? d[0] : d[0] + 1;
}

// Split a rectangle into 4 rectangles.
void decomp_rect_sizes (const Int nx, const Int ny, Int xd[2], Int yd[2]) {
  split(nx, xd);
  split(ny, yd);
  if (xd[1] > xd[0] && yd[1] > yd[0])
    std::swap(yd[0], yd[1]);
}

// Decompose a rectangle with coords into up to 4 non-empty pieces.
// p is (x y dx dy).
void decomp_rect (const Int p[4], Int kids[4][4], Int& nkids) {
  Int xd[2], yd[2];
  decomp_rect_sizes(p[2], p[3], xd, yd);
  nkids = 0;
  Int osx = p[0];
  for (Int i = 0; i < 2; ++i) {
    Int osy = p[1];
    for (Int j = 0; j < 2; ++j) {
      if (xd[i] > 0 && yd[j] > 0) {
        Int* kid = kids[nkids];
        kid[0] = osx;
        kid[1] = osy;
        kid[2] = xd[i];
        kid[3] = yd[j];
        ++nkids;
      }
      osy += yd[j];
    }
    osx += xd[i];
  }
}

// Compute the list of cells on a face covered by a rectangle.
void p2xy (const Int ne, const Int p[4], Int* xy, Int& nxy) {
  assert(p[2]*p[3] <= 4);
  nxy = 0;
  for (Int j = 0; j < p[3]; ++j) {
    const Int os = (p[1] + j)*ne;
    for (Int i = 0; i < p[2]; ++i)
      xy[nxy++] = os + p[0] + i;
  }
}

// Is this rectangle sufficiently small to form a leaf in the tree?
inline bool is_leaf (const Int p[4]) { return p[2] <= 2 && p[3] <= 2; }

// decomp_face helper.
Int decomp_face_recurse (const Int ne, const Int n_user_slots,
                         std::vector<Int>& tree, const Int p[4],
                         const Int parent) {
  const Int slot = tree.size();
  if (is_leaf(p)) {
    Int xy[4], nxy;
    p2xy(ne, p, xy, nxy);
    tree.resize(tree.size() + nxy + 2 + n_user_slots);
    tree[slot  ] = nxy;
    tree[slot+1] = parent;
    for (Int i = 0; i < nxy; ++i)
      tree[slot + i + 2] = xy[i];
  } else {
    Int kids[4][4], nkids;
    decomp_rect(p, kids, nkids);
    tree.resize(tree.size() + nkids + 2 + n_user_slots);
    tree[slot  ] = nkids;
    tree[slot+1] = parent;
    for (Int i = 0; i < nkids; ++i) {
      const Int kid_slot = decomp_face_recurse(ne, n_user_slots, tree, kids[i],
                                               slot);
      tree[slot + i + 2] = -kid_slot;
    }
  }
  return slot;
}

// Form a tree over the cells in a face of the cube.
void decomp_face (const Int ne, const Int n_user_slots,
                  std::vector<Int>& tree) {
  tree.clear();
  const Int root[] = {0, 0, ne, ne};
  decomp_face_recurse(ne, n_user_slots, tree, root, 0);
}

// Add index offsets to the tree's nodes and cells. This tree is over a face,
// and the offsets let multiple face trees to be put together to get a tree over
// the cube.
void add_offsets (Int* tree, const Int k,
                  const Int tos,    // Offset of node slot.
                  const Int cos)    // Offset of cell index.
{
  tree[k+1] += tos;
  for (Int i = 0, nkids = tree[k]; i < nkids; ++i) {
    Int& tki1 = tree[k+i+2];
    if (tki1 < 0) {
      add_offsets(tree, -tki1, tos, cos);
      tki1 -= tos;
    } else
      tki1 += cos;
  }
}

Int check_tree_recurse (const Int* tree, const Int k, const Int parent,
                        std::vector<Int>& cell_hits) {
  Int nerr = 0;
  if (tree::node_parent(tree, k) != parent)
    ++nerr;
  if (tree::node_nkids(tree, k) > 4)
    ++nerr;
  for (Int i = 0, nkids = tree::node_nkids(tree, k); i < nkids; ++i) {
    const Int e = tree::node_kid(tree, k, i);
    if (tree::node_has_cells(tree, k))
      ++cell_hits[e];
    else
      nerr += check_tree_recurse(tree, e, k, cell_hits);
  }
  return nerr;
}

// Unit test.
bool check_tree (const Int ncell, const Int* tree, const Int tree_size) {
  std::vector<Int> cell_hits(ncell, 0);
  const Int nerr = check_tree_recurse(tree, 0, 0, cell_hits);
  SIQK_STDERR_IF( nerr != 0, nerr << " ordinals didn't match.");
  bool all1 = true;
  for (Int i = 0; i < ncell; ++i)
    if (cell_hits[i] != 1) {
      all1 = false;
      break;
    }
  SIQK_STDERR_IF( ! all1, "Not all cells were covered, or a cell was covered"
                  << " more than once.");
  if ( ! all1) {
    siqk::prarr("tree", tree, tree_size);
    siqk::prarr("cell_hits", cell_hits.data(), cell_hits.size());
  }
  return all1;
}
} // namespace impl

bool make_cubedsphere_tree_over_cells (
  const Int ne, AIdxArray& tree, const Int n_user_slots,
  const bool perform_checks)
{
  std::vector<Int> face_tree;
  impl::decomp_face(ne, n_user_slots, face_tree);

  const Int
    face_tree_size = face_tree.size(),
    root_node_size = 5 + n_user_slots,
    next_node_size = 4 + n_user_slots,
    tree_top_size = root_node_size + 3*next_node_size;
  resize(tree, tree_top_size + 6*face_tree.size());
  // The root and next level's nodes are ((-y +x) (+y +z) (-x -z)), where each
  // entry is a face. Each face is then an index translation of face_tree.
  const Int face_os[] = {0, 1, 2, 4, 3, 5};
  
  // Root.
  tree[0] = 3;
  tree[1] = 0;
  for (Int i = 0; i < 3; ++i)
    tree[2+i] = -(root_node_size + i*next_node_size);
  // Root's children, each holding two adjacent faces, and such that all three
  // children are adjacent to each other.
  Int parents[6];
  for (Int i = 0; i < 3; ++i) {
    parents[i] = root_node_size + next_node_size*i;
    tree[root_node_size + next_node_size*i    ] =  2;
    tree[root_node_size + next_node_size*i + 1] =  0;
    tree[root_node_size + next_node_size*i + 2] =
      -(tree_top_size + face_os[2*i  ]*face_tree_size);
    tree[root_node_size + next_node_size*i + 3] =
      -(tree_top_size + face_os[2*i+1]*face_tree_size);
  }
  // Each face.
  for (Int i = 0; i < 6; ++i) {
    const Int tos = tree_top_size + face_os[i]*face_tree_size;
    const Int cos = face_os[i]*ne*ne;
    Int* const t = tree.data() + tos;
    memcpy(t, face_tree.data(), face_tree.size()*sizeof(Int));
    impl::add_offsets(t, 0, tos, cos);
    t[1] = parents[i/2];
  }

  if (perform_checks)
    return (impl::check_tree(ne*ne, face_tree.data(), face_tree.size()) &&
            impl::check_tree(6*ne*ne, tree.data(), nslices(tree)));

  return true;
}

namespace {
inline Int get_cube_face_idx (const Real& x, const Real& y, const Real& z) {
  const Real ax = std::abs(x), ay = std::abs(y), az = std::abs(z);
  if (ax >= ay) {
    if (ax >= az) return x > 0 ? 1 : 3;
    else return z > 0 ? 4 : 5;
  } else {
    if (ay >= az) return y > 0 ? 2 : 0;
    else return z > 0 ? 4 : 5;
  }
}

inline void map_sphere_coord_to_face_coord (
  const Int& face_idx, const Real& x, const Real& y, const Real& z,
  Real& fx, Real& fy)
{
  static constexpr Real theta_max = 0.25*M_PI;
  Real d;
  switch (face_idx) {
  case  0: d = std::abs(y); fx =  x/d; fy = z/d; break;
  case  1: d = std::abs(x); fx =  y/d; fy = z/d; break;
  case  2: d = std::abs(y); fx = -x/d; fy = z/d; break;
  case  3: d = std::abs(x); fx = -y/d; fy = z/d; break;
  case  4: d = std::abs(z); fx =  x/d; fy = y/d; break;
  default: d = std::abs(z); fx = -x/d; fy = y/d; break;
  }
  fx = std::atan(fx) / theta_max;
  fy = std::atan(fy) / theta_max;
}
} // namespace

Int get_cell_idx (const Int ne, const Real angle, const Real* const R,
                  Real x, Real y, Real z) {
  if (angle != 0) {
    // R'(x,y,z) to bring the point to the unrotated grid.
    Real v[3] = {x,y,z};
    x = R[0]*v[0] + R[3]*v[1] + R[6]*v[2];
    y = R[1]*v[0] + R[4]*v[1] + R[7]*v[2];
    z = R[2]*v[0] + R[5]*v[1] + R[8]*v[2];
  }
  const Int face_idx = get_cube_face_idx(x, y, z);
  Real fx, fy;
  map_sphere_coord_to_face_coord(face_idx, x, y, z, fx, fy);
  Int ix, iy;
  {
    iy = static_cast<Int>(std::floor(0.5*(1 + fy)*ne));
    if (iy < 0) iy = 0;
    else if (iy >= ne) iy = ne-1;
    ix = static_cast<Int>(std::floor(0.5*(1 + fx)*ne));
    if (ix < 0) ix = 0;
    else if (ix >= ne) ix = ne-1;
  }
  const Int face_cell_idx = (ne*iy + ix);
  assert(face_cell_idx < ne*ne);
  return ne*ne*face_idx + face_cell_idx;
}

void make_nonuniform (AVec3s& geo_p) {
  const Int n = nslices(geo_p);
  const Real axis[] = {1,1,1};
  const Real f = 0.5;
  Real R[9];
  form_rotation(axis, 0.2*M_PI, R);
  for (int i = 0; i < n; ++i) {
    auto p = slice(geo_p,i);
    Real p0[3];
    for (int d = 0; d < 3; ++d) p0[d] = p[d];
    for (int d = 0; d < 3; ++d) {
      Real a = 0;
      for (int j = 0; j < 3; ++j) a += R[3*j+d]*p0[j];
      p[d] = a;
    }
    for (int d = 1; d < 3; ++d) p[d] *= f;
    const Real r = std::sqrt(square(p[0]) + square(p[1]) + square(p[2]));
    for (int d = 0; d < 3; ++d) p[d] /= r;
    for (int d = 0; d < 3; ++d) p0[d] = p[d];
    for (int d = 0; d < 3; ++d) {
      Real a = 0;
      for (int j = 0; j < 3; ++j) a += R[3*d+j]*p0[j];
      p[d] = a;
    }
  }
}

void rotate_grid (const Real axis[3], const Real angle, Real* R, AVec3s& geo_p) {
  form_rotation(axis, angle, R);
  const Int n = nslices(geo_p);
  for (Int i = 0; i < n; ++i) {
    auto p = slice(geo_p,i);
    Real p0[3];
    for (int d = 0; d < 3; ++d) p0[d] = p[d];
    for (Int d = 0; d < 3; ++d) {
      const Real* const row = R + 3*d;
      p[d] = siqk::SphereGeometry::dot(row, p0);
    }
    siqk::SphereGeometry::normalize(p);
  }
}

} // namespace mesh
} // namespace slmm
