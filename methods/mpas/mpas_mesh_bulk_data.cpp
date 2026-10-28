#ifdef MPAS_HAVE_NETCDF
#include "mpas_mesh_bulk_data.hpp"

#include <netcdf>
using namespace netCDF;

/* Indexing conventions are a bit messy. The nc file and F90 use base-1 for
   everything. Thus, MpasMesh has base-1 indexing internally, but getter
   functions convert to base 0. In the following, c(l,g)id work arrays use a mix
   of base-0 and base-1, c ids are base-1 in MpasMesh, and e and v l,glids use
   base-1 for everything.
 */
struct MpasMeshBulkData::Impl {
  struct BulkData {
    std::vector<int> c2ne, c2c, c2e, c2v, e2c, e2v, v2c, v2e, e2e, e2ne;
    std::vector<Real> e2w, xc, yc, zc, xe, ye, ze, xv, yv, zv,
      dvEdge, dcEdge, angleEdge, areaCell;
  };

  bool ok;
  MpasMesh::Ptr m;
  BulkData bd;

  void read (const std::string& fname, Distribution& d, const Int nlevel) {
    ok = false;
    NcFile f(fname, NcFile::read);

    const auto dim = [&] (const char* name) { return f.getDim(name).getSize(); };
    const Int ncell = dim("nCells"), nedge = dim("nEdges"), nvtx = dim("nVertices"),
      maxedge = dim("maxEdges"), maxedge2 = dim("maxEdges2"), vtxdeg = dim("vertexDegree");
    const Int nlvl = nlevel < 0 ? dim("nVertLevels") : nlevel;

    const Int nlcell = d.c.gids.m();
    Int nledge, nlvtx;
    std::map<Int,Int> cg2lid, eg2lid, vg2lid;
    std::vector<Int> egids, vgids;
    Int neown, nvown;

    {
      std::vector<int> vi;
      const auto ivar = [&] (const char* name, const int& n) {
        if (static_cast<Int>(vi.size()) < n) vi.resize(n);
        f.getVar(name).getVar(vi.data());
      };

      ivar("nEdgesOnCell", ncell);
      bd.c2ne.resize(nlcell);
      for (Int i = 0; i < nlcell; ++i) bd.c2ne[i] = vi[d.c.gids[i]];

      // Order LIDs of x items so that any that is adjacent to an owned cell
      // comes first.
      const auto calc_c2x = [&] (std::vector<Int>& c2x, std::vector<Int>& xgids,
                                 std::map<Int,Int>& g2l, Int& nox) {
        c2x.resize(nlcell*maxedge);
        std::set<Int> uoxgids; // adjacent to owned cell
        for (Int i = 0; i < d.c.nown; ++i)
          for (Int j = 0; j < bd.c2ne[i]; ++j) {
            const Int x = vi[d.c.gids[i]*maxedge + j];
            uoxgids.insert(x);
          }
        nox = static_cast<Int>(uoxgids.size());
        xgids.reserve(nox);
        for (const auto& x : uoxgids) xgids.push_back(x);
        std::set<Int> uhxgids;
        for (Int i = d.c.nown; i < nlcell; ++i)
          for (Int j = 0; j < bd.c2ne[i]; ++j) {
            const Int x = vi[d.c.gids[i]*maxedge + j];
            if (uoxgids.find(x) == uoxgids.end()) uhxgids.insert(x);
          }
        xgids.reserve(xgids.size() + uhxgids.size());
        for (const auto& x : uhxgids) xgids.push_back(x);
        for (size_t i = 0; i < xgids.size(); ++i) g2l[xgids[i]] = i+1;
        for (Int i = 0; i < nlcell; ++i)
          for (Int j = 0; j < bd.c2ne[i]; ++j)
            c2x[i*maxedge + j] = g2l[vi[d.c.gids[i]*maxedge + j]];
      };
      ivar("edgesOnCell", ncell*maxedge);
      calc_c2x(bd.c2e, egids, eg2lid, neown);
      init_my_edge_data(maxedge, bd.c2e.data(), bd.c2ne.data(), egids.data(),
                        eg2lid, neown, d);
      // Could be done elsewhere, but for caller's convenience, do it in this
      // routine.
      exchange_edge_data(d);
      nledge = egids.size();
      ivar("verticesOnCell", ncell*maxedge);
      calc_c2x(bd.c2v, vgids, vg2lid, nvown);
      nlvtx = vgids.size();

      for (Int i = 0; i < nlcell; ++i) cg2lid[d.c.gids[i]+1] = i+1;

      const auto calc_e2x = [&] (const std::map<Int,Int>& g2l, std::vector<Int>& e2x) {
        e2x.resize(2*nledge);
        for (Int lei = 0; lei < nledge; ++lei) {
          const Int geim1 = egids[lei] - 1;
          auto it = g2l.find(vi[2*geim1]);
          e2x[2*lei  ] = it == g2l.end() ? -1 : it->second;
          it = g2l.find(vi[2*geim1+1]);
          e2x[2*lei+1] = it == g2l.end() ? -1 : it->second;
        }      
      };
      ivar("cellsOnEdge", 2*nedge); calc_e2x(cg2lid, bd.e2c);
      ivar("verticesOnEdge", 2*nedge); calc_e2x(vg2lid, bd.e2v);

      const auto calc_v2x = [&] (const std::map<Int,Int>& g2l, std::vector<Int>& v2x) {
        v2x.resize(vtxdeg*nlvtx);
        for (Int lvi = 0; lvi < nlvtx; ++lvi) {
          const Int gvim1 = vgids[lvi] - 1;
          for (Int j = 0; j < vtxdeg; ++j) {
            auto it = g2l.find(vi[vtxdeg*gvim1 + j]);
            v2x[vtxdeg*lvi + j] = it == g2l.end() ? -1 : it->second;
          }
        }
      };
      ivar("cellsOnVertex", nvtx*vtxdeg);
      calc_v2x(cg2lid, bd.v2c);
      ivar("edgesOnVertex", nvtx*vtxdeg);
      calc_v2x(eg2lid, bd.v2e);

      ivar("nEdgesOnEdge", nedge);
      bd.e2ne.resize(nledge);
      for (Int lei = 0; lei < nledge; ++lei) bd.e2ne[lei] = vi[egids[lei]-1];

      ivar("edgesOnEdge", nedge*maxedge2);
      bd.e2e.resize(nledge*maxedge2);
      for (Int lei = 0; lei < nledge; ++lei) {
        const Int geim1 = egids[lei] - 1;
        for (Int j = 0; j < bd.e2ne[lei]; ++j)
          bd.e2e[lei*maxedge2 + j] = eg2lid[vi[geim1*maxedge2 + j]];
      }

      ivar("cellsOnCell", ncell*maxedge);
      bd.c2c.resize(nlcell*maxedge);
      assert(d.c.gids.m() == nlcell);
      for (Int lci = 0; lci < nlcell; ++lci) {
        const Int gcim1 = d.c.gids[lci]; // base-0
        assert(gcim1 >= 0 && gcim1 < ncell);
        for (Int j = 0; j < bd.c2ne[lci]; ++j) {
          const auto it = cg2lid.find(vi[gcim1*maxedge + j]);
          const Int lcj = it == cg2lid.end() ? -1 : it->second;
          assert(lcj == -1 || (lcj > 0 && lcj <= nlcell)); // base-1
          bd.c2c[lci*maxedge + j] = lcj;
        }
      }
    }

    {
      std::vector<Real> vr;
      const auto rvar = [&] (const char* name, const int& n) {
        if (static_cast<Int>(vr.size()) < n) vr.resize(n);
        f.getVar(name).getVar(vr.data());
      };

      const auto rsub = [&] (const char* name, const Int gn, const std::vector<Int>& gids,
                             std::vector<Real>& v) {
        rvar(name, gn);
        v.reserve(gids.size());
        for (const auto& gid : gids) v.push_back(vr[gid-1]);
      };

      // Base-1 indexing.
      std::vector<Int> cgids(d.c.gids.m());
      for (Int i = 0; i < d.c.gids.m(); ++i) cgids[i] = d.c.gids[i] + 1;

      rsub("xCell", ncell, cgids, bd.xc);
      rsub("yCell", ncell, cgids, bd.yc);
      rsub("zCell", ncell, cgids, bd.zc);
      rsub("xEdge", nedge, egids, bd.xe);
      rsub("yEdge", nedge, egids, bd.ye);
      rsub("zEdge", nedge, egids, bd.ze);
      rsub("xVertex", nvtx, vgids, bd.xv);
      rsub("yVertex", nvtx, vgids, bd.yv);
      rsub("zVertex", nvtx, vgids, bd.zv);
      rsub("dvEdge", nedge, egids, bd.dvEdge);
      rsub("dcEdge", nedge, egids, bd.dcEdge);
      rsub("angleEdge", nedge, egids, bd.angleEdge);
      rsub("areaCell", ncell, cgids, bd.areaCell);

      rvar("weightsOnEdge", nedge*maxedge2);
      bd.e2w.resize(nledge*maxedge2);
      for (Int lei = 0; lei < nledge; ++lei) {
        const Int geim1 = egids[lei] - 1;
        for (Int j = 0; j < maxedge2; ++j)
          bd.e2w[lei*maxedge2 + j] = vr[geim1*maxedge2 + j];
      }
    }
    
#define dat(a) bd.a.data()
    m = std::make_shared<MpasMesh>(
      nlcell, nledge, nlvtx, maxedge, vtxdeg,
      dat(c2ne), dat(c2c), dat(c2e), dat(c2v), dat(e2c), dat(e2v), dat(v2c), dat(v2e),
      dat(xc), dat(yc), dat(zc), dat(xe), dat(ye), dat(ze), dat(xv), dat(yv), dat(zv),
      dat(dvEdge), dat(dcEdge), dat(angleEdge), dat(areaCell),
      dat(e2e), dat(e2ne), dat(e2w), maxedge2,
      nlvl);
#undef dat
    init_extras(*m);
    m->ncown = d.c.nown; m->neown = neown; m->nvown = nvown;
    ok = true;
  }

  void check (const Distribution& d) {
    std::stringstream ss;
    ss << "MpasMeshBulkData::check (";
    ss << d.p->rank() << "):";
    const std::string slbl = ss.str();
    const char* lbl = slbl.c_str();
    if ( ! ok) { printf("%s !ok on entry\n", lbl); return; }
    ok = true;
    const auto& t = m->t;
    const auto& dist = d.cdhalo;
    for (Int ci = 0; ci < t.ncell; ++ci)
      for (Int i = 0; i < t.c2ne[ci]; ++i) {
        const auto ei = t.get_c2e(ci, i);
        bool fnd = false;
        for (Int j = 0; j < 2; ++j)
          if (t.get_e2c(ei, j) == ci) fnd = true;
        if ( ! fnd) ok = false;
      }
    if ( ! ok) { printf("%s !ok c2e, e2c\n", lbl); return; }
    for (Int ci = 0; ci < t.ncell; ++ci) {
      if (dist[ci] == 1) continue;
      for (Int i = 0; i < t.c2ne[ci]; ++i) {
        const auto ei = t.get_c2e(ci, i);
        if (ei < 0) ok = false;
        for (Int j = 0; j < 2; ++j)
          if (t.get_e2c(ei, j) < 0) ok = false;
        const auto vi = t.get_c2v(ci, i);
        for (Int j = 0; j < t.vdegree; ++j) {
          if (t.get_v2c(vi, j) < 0) ok = false;
          if (t.get_v2e(vi, j) < 0) ok = false;
        }
        for (Int j = 0; j < t.e2ne[ei]; ++j)
          if (t.get_e2e(ei, j) < 0) ok = false;
      }
    }
    if ( ! ok) { printf("%s !ok dist != 1\n", lbl); return; }
    for (Int ci = 0; ci < t.ncell; ++ci)
      for (Int i = 0; i < t.c2ne[ci]; ++i) {
        if (t.get_c2v(ci, i) < 0) ok = false;
        if (t.get_c2e(ci, i) < 0) ok = false;
      }
    if ( ! ok) { printf("%s !ok c2v\n", lbl); return; }
    for (Int ci = 0; ci < t.ncell; ++ci) {
      Int n = 0;
      for (Int j = 0; j < t.c2ne[ci]; ++j) {
        const Int cj = t.get_c2c(ci, j);
        if (cj < 0) continue;
        ++n;
        bool fnd = false;
        for (Int k = 0; k < t.c2ne[cj]; ++k)
          if (t.get_c2c(cj, k) == ci) fnd = true;
        if ( ! fnd) ok = false;
      }
      if (n < 2) ok = false;
    }
    if ( ! ok) { printf("%s !ok c2c\n", lbl); return; }
  }
};

MpasMeshBulkData::MpasMeshBulkData () {
  p = std::make_shared<Impl>();
  p->ok = false;
}

bool MpasMeshBulkData::is_ok () const { return p->ok; }

MpasMesh::Ptr MpasMeshBulkData::get_mpas_mesh () { return p->m; }

MpasMeshBulkData::Ptr MpasMeshBulkData
::create_and_update_distribution (const std::string& fname, Distribution& d,
                                  const Int nlevel) {
  const auto nm = std::make_shared<MpasMeshBulkData>();
  nm->p->read(fname, d, nlevel);
#ifndef NDEBUG
  nm->p->check(d);
#endif
  if ( ! nm->is_ok()) return nullptr;
  return nm;
}

#endif // MPAS_HAVE_NETCDF
