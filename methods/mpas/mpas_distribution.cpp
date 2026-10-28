#include "mpas_distribution.hpp"
#include "mpas_util.hpp"

#include <set>
#include <map>

// Cell l,glids are base-0 in this file.

static void
collect_own_and_halo_cells (const std::vector<Int>& c2r, const CellGraph& cg,
                            const Int mypid, const Int szhalo,
                            std::vector<char>& in /* dhalo */,
                            std::vector<Int>& v /* cgid_own_halo */) {
  const Int ncell = c2r.size();
  assert(cg.ncell() == ncell);
  const Int notin = -1;
  in.clear();
  in.resize(ncell, notin);
  // Own cells.
  Int psym = 0;
  for (Int ci = 0; ci < ncell; ++ci)
    if (c2r[ci] == mypid)
      in[ci] = psym;
  // Breadth-first search.
  for (Int halo = 0; halo < szhalo; ++halo) {
    const Int nsym = szhalo - halo;
    for (Int ci = 0; ci < ncell; ++ci) {
      // Check only the current frontier cells.
      if (in[ci] != psym) continue;
      // Add any un-included nbrs of this frontier cell. Mark it as a new
      // frontier cell.
      for (Int j = cg.cells[ci]; j < cg.cells[ci+1]; ++j) {
        const Int cin = cg.nbrs[j];
        if (in[cin] == -1) in[cin] = nsym;
      }
    }
    psym = nsym;
  }
  // Get number of owned + halo cells.
  Int n = 0;
  for (Int ci = 0; ci < ncell; ++ci)
    if (in[ci] != notin)
      ++n;
  // Make a list of these cells.
  v.clear();
  v.reserve(n);
  for (Int ci = 0; ci < ncell; ++ci)
    if (in[ci] != notin)
      v.push_back(ci);
}

static void
exchange_lids (const mpi::Parallel& p, const Aci1& pids, const Int nown,
               const Int nhalo, const Aci1& pidptrs, const Aci1& gids,
               Ai1& rlids, Ai1& olids_ptrs, Ai1& olids) {
  if (p.size() == 1) {
    assert(nhalo == 0);
    olids_ptrs.reset(2);
    olids_ptrs[0] = olids_ptrs[1] = 0;
    return;
  }
  const Int tag = 109;
  const Int npid = pids.m();
  std::vector<mpi::Request> rreq(npid), sreq(npid);
  auto& rcnt = olids_ptrs;
  rcnt.reset(npid+1);
  { // Exchange counts.
    std::vector<Int> isbuf(npid);  
    for (Int pi = 1; pi < npid; ++pi) {
      const Int pid = pids(pi);
      // Set up recv of remote #lids needed.
      mpi::irecv(p, &rcnt[pi], 1, pid, tag, &rreq[pi]);
      // Send remotes the #lids I need.
      isbuf[pi] = pidptrs(pi+1) - pidptrs(pi);
      mpi::isend(p, &isbuf[pi], 1, pid, tag, &sreq[pi]);
    }
    mpi::waitall(npid-1, &rreq[1]);
    mpi::waitall(npid-1, &sreq[1]);
  }
  // Allocate an extra buffer using counts.
  Int ntotcnt = 0;
  for (Int pi = 1; pi < npid; ++pi) ntotcnt += rcnt[pi];
  auto& ibuf = olids;
  ibuf.reset(ntotcnt);
  ibuf.set(-1);
  // Send and recv gids requests.
  for (Int pi = 1, i0 = 0; pi < npid; ++pi) {
    const Int pid = pids(pi);
    // Set up recv of remote gids requests.
    mpi::irecv(p, &ibuf[i0], rcnt[pi], pid, tag, &rreq[pi]);
    i0 += rcnt[pi];
    // Send gids requests.
    mpi::isend(p, &gids(pidptrs(pi)), pidptrs(pi+1) - pidptrs(pi),
               pid, tag, &sreq[pi]);
  }
  mpi::waitall(npid-1, &rreq[1]);
  mpi::waitall(npid-1, &sreq[1]);
  // Map gid to my lid.
  std::map<Int,Int> g2l;
  for (Int i = 0; i < nown; ++i) g2l[gids[i]] = i;
  for (Int i = 0; i < ntotcnt; ++i) {
    assert(ibuf[i] >= 0);
    const auto it = g2l.find(ibuf[i]);
    assert(it != g2l.end());
    ibuf[i] = it->second;
  }
  // Send and recv lids to fill rlids.
  rlids.reset(nhalo);
#ifndef NDEBUG
  rlids.set(-1);
#endif
  rcnt[0] = 0;
  for (Int pi = 1, i0 = 0; pi < npid; ++pi) {
    const Int pid = pids(pi);
    // Set up recv of remote lids responses.
    mpi::irecv(p, &rlids(pidptrs(pi) - nown), pidptrs(pi+1) - pidptrs(pi),
               pid, tag, &rreq[pi]);
    // Send lids responses.
    mpi::isend(p, &ibuf[i0], rcnt[pi], pid, tag, &sreq[pi]);
    i0 += rcnt[pi];
  }
  // Shift rcnt entries and cumsum to make olids_ptrs.
  for (Int pi = npid-1; pi >= 0; --pi) rcnt[pi+1] = rcnt[pi];
  for (Int pi = 1; pi <= npid; ++pi) rcnt[pi] += rcnt[pi-1];
  assert(rcnt[npid] == ntotcnt);
  mpi::waitall(npid-1, &rreq[1]);
  mpi::waitall(npid-1, &sreq[1]);
#ifndef NDEBUG
  for (Int i = 0, n = rlids.m(); i < n; ++i)
    assert(rlids(i) >= 0);
  for (Int i = 0, n = olids.m(); i < n; ++i)
    assert(olids(i) >= 0 && olids(i) < nown);
#endif
}

static void exchange_clids (Distribution& d) {
  auto& c = d.c;
  exchange_lids(*d.p, d.pids, c.nown, c.nhalo, c.pidptrs, c.gids, c.rlids,
                c.olids_ptrs, c.olids);
}

Distribution::Ptr Distribution
::create (const mpi::Parallel::Ptr& p, const std::vector<Int>& c2r, const CellGraph& cg,
          const Int szhalo) {
  // since we use char for the dist array
  if (szhalo < 1 || szhalo >= 128) {
    if (p->amroot()) printf("Distribution::create: szhalo must be > 1 and <= 127\n");
    return nullptr;
  }

  auto d = std::make_shared<Distribution>();
  d->p = p;
  d->szhalo = szhalo;
  const Int mypid = p->rank();

  std::vector<Int> cgid_own_halo;
  std::vector<char> gcdhalo;
  collect_own_and_halo_cells(c2r, cg, mypid, szhalo, gcdhalo, cgid_own_halo);
  const Int ncell = cgid_own_halo.size();
  {
    // Count own and halo cells.
    Int nown = 0;
    std::set<Int> pids;
    for (Int i = 0; i < ncell; ++i) {
      const Int pid = c2r[cgid_own_halo[i]];
      pids.insert(pid);
      if (pid == mypid) ++nown;
    }
    d->c.nown = nown;
    d->c.nhalo = ncell - nown;
    // Collect remotes with whom I may interact.
    const Int nremote = pids.size();
    d->pids.reset(std::max(1, nremote));
    d->pids(0) = mypid;
    Int i = 1;
    for (const auto& r : pids)
      if (r != mypid)
        d->pids(i++) = r;
  }
  // Group cgids by pid, starting with me.
  const Int nremote = d->pids.m();
  d->c.gids.reset(ncell);
  d->c.pidptrs.reset(nremote+1);
  d->c.pidptrs(0) = 0;
  for (Int i = 0; i < nremote; ++i) {
    d->c.pidptrs(i+1) = d->c.pidptrs(i);
    const Int pid = d->pids(i);
    for (Int j = 0; j < ncell; ++j) {
      const Int ci = cgid_own_halo[j];
      if (c2r[ci] == pid)
        d->c.gids(d->c.pidptrs(i+1)++) = ci;
    }
  }
  assert(d->c.pidptrs(nremote) == ncell);
  // Get remote LIDs for halo cells.
  exchange_clids(*d);
  // Get distance to halo edge in cgcids order.
  d->cdhalo.reset(ncell);
  for (Int i = 0; i < ncell; ++i) d->cdhalo(i) = gcdhalo[d->c.gids(i)];

  return d;
}

/* Determine edge LIDs so that edges are grouped by owning proc. An edge is
   owned by two procs; if a proc does not own it but shares it, the proc chooses
   one source. Procs don't have to agree on which owns the edge; the only
   requirement is that a proc knows where to get an unowned edge's data if it
   doesn't have it already. Thus, neown summed across procs is expected to be >
   than the number of edges.
*/
void init_my_edge_data (const Int maxedge, Int* const c2e, const Int* c2ne,
                        Int* egids, std::map<Int,Int>& eg2l, Int& neown,
                        Distribution& d) {
  const Int nremote = d.pids.m(), ncell = d.c.nown + d.c.nhalo, nedge = eg2l.size();
  auto& e = d.e;

  // Map cell LID-1 to remote LID.
  std::vector<Int> c2r(ncell, -1);
  for (Int r = 0; r < nremote; ++r)
    for (Int c = d.c.pidptrs[r]; c < d.c.pidptrs[r+1]; ++c)
      c2r[c] = r;

  // Map prelim edge LID to cell (LID-1)s. 'pe' means 'preliminary edge LID',
  // the edge LIDs on input.
  std::vector<std::array<Int,2> > pe2cs(nedge+1); // slot 0 is unused
  for (Int i = 1; i <= nedge; ++i) pe2cs[i][0] = pe2cs[i][1] = -1;
  for (Int c = 0; c < ncell; ++c)
    for (Int j = 0; j < c2ne[c]; ++j) {
      const Int lpe = c2e[c*maxedge + j];
      assert(lpe < 0 || lpe > 0);
      if (lpe < 0) break;
      assert(lpe <= nedge);
      if (pe2cs[lpe][0] < 0)
        pe2cs[lpe][0] = c;
      else {
        assert(pe2cs[lpe][1] == -1);
        pe2cs[lpe][1] = c;
      }
    }

  // Map prelim edge LID to remote LID.
  std::vector<Int> pe2r(nedge+1, -1); // slot 0 is unused
  Int slot = 0;
  for (Int c = 0; c < ncell; ++c)
    for (Int j = 0; j < c2ne[c]; ++j) {
      const Int lpe = c2e[c*maxedge + j];
      if (lpe < 0) break;
      if (pe2r[lpe] >= 0) continue;
      if (c < d.c.nown) {
        assert(lpe <= neown);
        const Int r = c2r[c];
        assert(r == 0);
        pe2r[lpe] = r;
      } else if (pe2cs[lpe][slot] < 0) {
        const Int c = pe2cs[lpe][(slot+1)%2];
        assert(c >= 0 && c < ncell);
        const Int r = c2r[c];
        assert(r >= 0);
        pe2r[lpe] = r;
      } else {
        // Arbitrarily choose owning remote.
        const Int c = pe2cs[lpe][slot];
        assert(c >= 0 && c < ncell);
        const Int r = c2r[c];
        pe2r[lpe] = r;
        slot = (slot + 1) % 2;
      }
    }

  const auto init_epidptrs = [&] () {
    for (Int i = 0; i <= nremote; ++i) e.pidptrs[i] = 0;
    for (Int i = 1; i <= nedge; ++i) {
      assert(pe2r[i] >= 0);
      e.pidptrs[pe2r[i]+1]++;
    }
    for (Int i = 1; i <= nremote; ++i) e.pidptrs[i] += e.pidptrs[i-1];
  };
  
  // Determine new edge LIDs, grouped by remotes.
  e.pidptrs.reset(nremote+1);
  e.gids.reset(nedge);
  for (Int i = 0; i < nedge; ++i) e.gids[i] = -1;
  init_epidptrs();
  for (Int i = 1; i <= nedge; ++i) {
    const Int egid = egids[i-1];
    assert(egid >= 1);
    e.gids[e.pidptrs[pe2r[i]]++] = egid;
  }
  assert(e.pidptrs[nremote] == nedge);
  assert(e.pidptrs[0] == neown);
  init_epidptrs();
  // Reset inputs using new edge LIDs.
  //   eg2l
  for (Int i = 0; i < nedge; ++i) {
    assert(e.gids[i] >= 1);
    eg2l[e.gids[i]] = i+1;
  }
  //   c2e
  for (Int c = 0; c < ncell; ++c)
    for (Int j = 0; j < c2ne[c]; ++j) {
      if (c2e[c*maxedge + j] < 0) break;
      const Int elid_prelim = c2e[c*maxedge + j];
      assert(elid_prelim > 0);
      const Int elid_new = eg2l[egids[elid_prelim-1]];
      assert(elid_new > 0);
      c2e[c*maxedge + j] = elid_new;
    }
  //   egids
  for (Int i = 0; i < nedge; ++i) {
    assert(e.gids[i] >= 1);
    egids[i] = e.gids[i];
  }
  e.nown = neown;
  e.nhalo = nedge - neown;
}

// Init relids, oelids_ptrs, oelids.
void exchange_edge_data (Distribution& d) {
  auto& e = d.e;
  exchange_lids(*d.p, d.pids, e.nown, e.nhalo, e.pidptrs, e.gids, e.rlids,
                e.olids_ptrs, e.olids);
}
