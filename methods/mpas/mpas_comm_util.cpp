#include "mpas_comm_util.hpp"
#include "mpas_geometry.hpp"

Int get_num_own_cells (const Topo& t, const Comm& he) {
  if ( ! he.parallel()) return t.ncell;
  const auto& dist = he.cdhalo();
  for (Int i = 0; i < dist.m(); ++i)
    if (dist(i) != 0)
      return i;
  assert(0);
  return -1;
}

Real get_avg_cell_spacing (const MpasMesh& m) {
  Real s = 0;
  for (int i = 0; i < m.t.nedge; ++i) s += m.g.dcEdge[i];
  return s / m.t.nedge;
}

Real get_avg_radius (const MpasMesh& m) {
  Real R = 0;
  for (int i = 0; i < m.t.ncell; ++i) {
    Vec cc;
    fill(m.c, 0, cc);
    R += geo::norm(cc);
  }
  return R/m.t.ncell;
}

void get_avg_radius_and_cell_spacing (const MpasMesh& m, Comm& he,
                                      Real& radius, Real& cell_spacing) {
  Real s = 0, R = 0;
  for (int i = 0; i < m.t.nedge; ++i)
    s += m.g.dcEdge[i];
  for (int i = 0; i < m.t.ncell; ++i) {
    Vec cc;
    fill(m.c, i, cc);
    R += geo::norm(cc);
  }
  Real laccum[4] = {s, Real(m.t.nedge), R, Real(m.t.ncell)}, gaccum[4];
  he.all_reduce(Comm::SUM, sizeof(laccum)/sizeof(*laccum), laccum, gaccum);
  cell_spacing = gaccum[0]/gaccum[1];
  radius = gaccum[2]/gaccum[3];
}

#ifndef MPAS_HAVE_NETCDF
CellGraph::ConstPtr create_CellGraph (const std::string& mesh_fname) {
  printf("CellGraph::ConstPtr create from netcdf unsupported b/c "
         "ndef MPAS_HAVE_NETCDF.\n");
  return nullptr;
}
#else
#include <netcdf>
using namespace netCDF;

CellGraph::ConstPtr create_CellGraph (const std::string& mesh_fname) {
  if ( ! fileexists(mesh_fname)) return nullptr;
  NcFile f(mesh_fname, NcFile::read);
  const auto dim = [&] (const char* name) { return f.getDim(name).getSize(); };
  const Int ncell = dim("nCells"), maxedges = dim("maxEdges");
  std::vector<int> c2c(ncell*maxedges);
  f.getVar("cellsOnCell").getVar(c2c.data());
  std::vector<int> c2ne(ncell);
  f.getVar("nEdgesOnCell").getVar(c2ne.data());
  const auto cg = std::make_shared<CellGraph>();
  const auto init_cells = [&] () {
    cg->cells.resize(ncell+1);
    cg->cells[0] = 0;
    for (Int ci = 0; ci < ncell; ++ci)
      cg->cells[ci+1] = cg->cells[ci] + c2ne[ci];
  };
  init_cells();
  cg->nbrs.resize(cg->cells[ncell]);
  for (Int ci = 0; ci < ncell; ++ci)
    for (Int j = 0; j < c2ne[ci]; ++j) {
      const Int cj = c2c[ci*maxedges + j] - 1; // to base-0
      assert(cj >= 0);
      assert(cg->cells[ci] < cg->cells[ncell]);
      cg->nbrs[cg->cells[ci]++] = cj;
    }
  init_cells();
  return cg;
}
#endif
