#ifndef INCLUDE_MPAS_MESH_HPP
#define INCLUDE_MPAS_MESH_HPP

#include "mpas_types.hpp"
#include "mpas_siqp.hpp"

struct ZetaRep;
struct Comm;

static constexpr int max_vtx = 9;
static constexpr int max_vdegree = 3;

static constexpr int max_1halo_slots = max_vtx+1;
static constexpr int max_2halo_slots = 24;
static constexpr int max_3halo_slots = 40;
static constexpr int max_4halo_slots = 64;
static constexpr int max_5halo_slots = 96;
static constexpr int max_6halo_slots = 147;
static constexpr int max_halo_slots[] = {max_1halo_slots, max_2halo_slots, max_3halo_slots,
                                         max_4halo_slots, max_5halo_slots, max_6halo_slots};
static constexpr int max_nhalo = sizeof(max_halo_slots)/sizeof(*max_halo_slots);

struct Topo {
  int ncell, nedge, nvtx, max_vtx, vdegree, max_e2w;
  const int* c2ne, * c2c, * c2e, * c2v, * e2c, * e2v, * v2c, * v2e, * e2e, * e2ne;
  const Real* e2w;

  Topo (int ncell_, int nedge_, int nvtx_, int max_vtx_, int vdegree_,
        const int* c2ne_, const int* c2c_, const int* c2e_, const int* c2v_,
        const int* e2c_, const int* e2v_,
        const int* v2c_, const int* v2e_,
        const int* e2e_, const int* e2ne_, const Real* e2w_, int max_e2w_)
    : ncell(ncell_), nedge(nedge_), nvtx(nvtx_), max_vtx(max_vtx_),
      vdegree(vdegree_), max_e2w(max_e2w_),
      c2ne(c2ne_), c2c(c2c_), c2e(c2e_), c2v(c2v_),
      e2c(e2c_), e2v(e2v_), v2c(v2c_), v2e(v2e_),
      e2e(e2e_), e2ne(e2ne_), e2w(e2w_)
  {}

  int get_c2c(int ci, int i) const {
    assert(ci >= 0 && ci < ncell);
    assert( i >= 0 &&  i < max_vtx);
    return c2c[max_vtx*ci + i] - 1;
  }

  int get_c2e(int ci, int ei) const {
    assert(ci >= 0 && ci < ncell);
    assert(ei >= 0 && ei < max_vtx);
    return c2e[max_vtx*ci + ei] - 1;
  }

  int get_c2v(int ci, int vi) const {
    assert(ci >= 0 && ci < ncell);
    // max_vtx is used for c2v, as well.
    assert(vi >= 0 && vi < max_vtx);
    return c2v[max_vtx*ci + vi] - 1;
  }

  int get_e2c(int ei, int ci) const {
    assert(ei >= 0 && ei < nedge);
    assert(ci >= 0 && ci < 2);
    return e2c[2*ei + ci] - 1;
  }

  int get_e2v(int ei, int vi) const {
    assert(ei >= 0 && ei < nedge);
    assert(vi >= 0 && vi < 2);
    return e2v[2*ei + vi] - 1;
  }

  int get_v2c(int vi, int ci) const {
    assert(vi >= 0 && vi < nvtx);
    assert(ci >= 0 && ci < vdegree);
    return v2c[vdegree*vi + ci] - 1;
  }

  int get_v2e(int vi, int ei) const {
    assert(vi >= 0 && vi < nvtx);
    assert(ei >= 0 && ei < vdegree);
    return v2e[vdegree*vi + ei] - 1;
  }

  int get_e2e(int ei, int i) const {
    assert(ei >= 0 && ei < nedge);
    assert(i >= 0 && i < max_e2w);
    return e2e[max_e2w*ei + i] - 1;
  }

  Real get_e2w(int ei, int i) const {
    assert(ei >= 0 && ei < nedge);
    assert(i >= 0 && i < max_e2w);
    return e2w[max_e2w*ei + i];
  }
};

struct Points {
  int n;
  const Real* x, * y, * z;
  Points (int n_, const Real* x_, const Real* y_, const Real* z_)
    : n(n_), x(x_), y(y_), z(z_)
  {}
};

inline void fill (const Points& p, const int i, Real v[3], bool normalize = false) {
  v[0] = p.x[i]; v[1] = p.y[i]; v[2] = p.z[i];
  if (normalize) {
    const auto len = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    v[0] /= len; v[1] /= len; v[2] /= len;
  }
}

struct CellCtrs : public Points {
  CellCtrs (int n, const Real* x, const Real* y, const Real* z) : Points(n, x, y, z) {}
};

struct EdgeCtrs : public Points {
  EdgeCtrs (int n, const Real* x, const Real* y, const Real* z) : Points(n, x, y, z) {}
};

struct Vtxs : public Points {
  Vtxs (int n, const Real* x, const Real* y, const Real* z) : Points(n, x, y, z) {}
};

struct MpasGeometry {
  const Real* dvEdge, * dcEdge, * angleEdge, * areaCell;
  std::vector<Real> my_angleEdge, my_areaCell;
  Array2D<Real> my_cellCentroid;
  std::vector<short> my_sign;

  MpasGeometry (const Real* dvEdge_, const Real* dcEdge_, const Real* angleEdge_,
                const Real* areaCell_)
    : dvEdge(dvEdge_), dcEdge(dcEdge_), angleEdge(angleEdge_), areaCell(areaCell_)
  {}
};

/* MPAS-O conventions:
     1. k=1 at ocean top.
     2. w = vertAleTransportTop is positive pointing to ocean surface.
     3. zMid(k+1) < zMid(k), consistent with w.
   SL internal conventions are the opposite of MPAS-O's for 2 and 3.
 */
struct Levels {
  int nmax;
  // (min,max)LevelCell and (min,max)LevelEdge(Top,Bot), where the index bottom
  // is the physical ocean top since level indices go from ocean top to bottom.
  const int* cmin, * cmax, * eminbot, * emaxbot, * emintop, * emaxtop;

  Levels(const int ncell, const int nedge, const int nlevel);
  // Set all cell-center zeta interface values to zint[0:nmax].
  void set_interface_z_all_cols(const Real* zint);
  // Set cell-center zeta interface vales to smoothly varying values. Each
  // column extends from zlo to zhi.
  void set_interface_spatially_varying_zeta_eul(const CellCtrs& c, const Real zlo,
                                                const Real zhi);

  // base-1 -> base-0 getters
  int get_cmin (int c) const { return cmin[c]-1; }
  int get_cmax (int c) const { return cmax[c]-1; }
  int get_eminbot (int e) const { return eminbot[e] - 1; }
  int get_emaxbot (int e) const { return emaxbot[e] - 1; }
  int get_emintop (int e) const { return emintop[e] - 1; }
  int get_emaxtop (int e) const { return emaxtop[e] - 1; }
  const Real* get_interface_z(int c) const;

  const Acr2& get_zint_array () const { return zint; }

private:
  // For standalone.
  std::vector<int> dmin, dmax;
  Acr2 zint;
};

struct MpasMesh {
  typedef std::shared_ptr<MpasMesh> Ptr;
  typedef std::shared_ptr<const MpasMesh> ConstPtr;

  Topo t;
  CellCtrs c;
  EdgeCtrs e;
  Vtxs v;
  MpasGeometry g;
  Levels l;
  std::shared_ptr<const ZetaRep> zeta_rep;

  // In the distributed case, ncown is the number of owned cells, which is <
  // t.ncell. neown is the number of edges, all ordered first, that are adjacent
  // to at least one owned cell, and similarly for vertices and nvown.
  Int ncown, neown, nvown;
  
  MpasMesh(int ncell, int nedge, int nvtx, int max_c2e, int vdegree,
           const PyInt* c2ne, const PyInt* c2c, const PyInt* c2e, const PyInt* c2v,
           const PyInt* e2c, const PyInt* e2v,
           const PyInt* v2c, const PyInt* v2e,
           const Real* xc, const Real* yc, const Real* zc,
           const Real* xe, const Real* ye, const Real* ze,
           const Real* xv, const Real* yv, const Real* zv,
           const Real* dvEdge, const Real* dcEdge, const Real* angleEdge, const Real* areaCell,
           // For edge velocity reconstruction
           const PyInt* e2e, const PyInt* e2ne, const Real* e2w, int max_e2w,
           const int nlevel = 1);
};

struct Stats {
  Real min, max, mean;
  int cnt;
  Stats () : min(1e30), max(-1e30), mean(0), cnt(0) {}
  void operator() (const Real& a) {
    min = std::min(min, a); max = std::max(max, a);
    mean += a;
    ++cnt;
  }
  void print (const std::string& label) {
    printf("%20s: %1.3e %1.3e %1.3e\n", label.c_str(), min, mean/cnt, max);
  }
};

siqp::sh::Mesh::Ptr create_siqp_mesh(const MpasMesh& mm);

// Init my_angleEdge, my_sign.
void init_extras(MpasMesh& m);
// Init my_areaCell, my_cellCentroid from m.
void init_more_extras(MpasMesh& mm, siqp::sh::Mesh& m);
// Init zeta_rep using Levels::get_interface_z.
void init_zeta_rep(const MpasMesh::Ptr& mm, const siqp::sh::Mesh& m,
                   const std::shared_ptr<Comm>& he);

inline int get_other_cell (const Topo& t, const int cme, const int ei) {
  const int ci0 = t.get_e2c(ei,0), ci1 = t.get_e2c(ei,1);
  assert(cme == ci0 || cme == ci1);
  return cme == ci0 ? ci1 : ci0;
}

void get_oriented_cell(const Topo& t, const Vtxs& v, const int ci,
                       int* vtxs, int* edges = nullptr);

void calc_tangent_plane_axes(Vec ctr, Vec xaxis, Vec yaxis);

void calc_xy(const Vec ctr, const Vec p, const Vec xaxis, const Vec yaxis,
             const bool scale_by_arc_length, Real& x, Real& y);

int map_cells_to_cell_plane(const MpasMesh& m, const int ci,
                            Vec xaxis, Vec yaxis, Real* xs, Real* ys,
                            int* stencil = nullptr,
                            const bool scale_by_arc_length = false,
                            const bool two_halo = false);
void map_edges_to_vertex_plane(const MpasMesh& m, const int vi,
                               Vec xaxis, Vec yaxis, Real* xs, Real* ys,
                               int* stencil = nullptr, const int nhop = 1,
                               const bool scale_by_arc_length = false);
void map_vtxs_to_cell_plane(const MpasMesh& m, const int ci,
                            Vec xaxis, Vec yaxis, Real* xs, Real* ys,
                            int* stencil = nullptr,
                            const bool scale_by_arc_length = false,
                            // In the case of a moving grid, replace m.c and m.v
                            // with the following cell centroid and vertices
                            // dimensioned as v_replacement(3,nlvl,nvtx).
                            const Real* cc_replacement = nullptr,
                            const Array3D<const Real>* v_replacement = nullptr,
                            const int v_lvl = -1);

void calc_cell_centroid(const MpasMesh& m, const int ci, Vec cc);

int augment_stencil(const MpasMesh& m, const int ci,
                    // ncurrent is the number of cells currently in stencil. The
                    // stencil array excludes ci.
                    int ncurrent,
                    // Add to the existing stencil. Return the total count.
                    int* stencil,
                    // n calls of augment_stencil is equivalent in output to
                    // setting nhalo = n.
                    const int nhalo = 1,
                    // In distributed case, cid can be -1 for a cell on the edge
                    // of the halo, so handle that case.
                    const bool permit_neg_cid = false);

#endif
