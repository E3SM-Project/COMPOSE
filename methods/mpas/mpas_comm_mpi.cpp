#include "mpas_comm.hpp"
#include "mpas_distribution.hpp"
#include "mpas_util.hpp"
#include "mpas_accum.hpp"
#include "mpas_mesh.hpp"

#include <cstdint>

// Hold an entity (LID, level) pair, such as cell LID and level. To communicate
// an LidLvlPair, reinterpret_cast it as an LidLvlPair::Id.
struct LidLvlPair {
  typedef uint32_t Id;
  static constexpr int sz = 10;
  LidLvlPair () {}
  LidLvlPair (Int lid, Int lvl) { set(lid, lvl); }
  void set (Int lid, Int lvl) { encode(lid, lvl, id_); }
  Int lid () const { return (Int) (id_ >> sz); }
  Int lvl () const { return (id_ << (32 - sz)) >> (32 - sz); }
  Id id () const { return id_; }
  static bool test();
  static void assert_sizes(int nent, int nlvl);
private:
  Id id_;
  static inline void encode(const Int lid, const Int lvl, Id& id);
};

// (int,int) -> Id
void LidLvlPair::encode (const Int lid, const Int lvl, uint& id) {
  assert(lid >= 0 && lvl >= 0);
  id = (((uint) lid) << sz) | ((uint) lvl);
}

bool LidLvlPair::test () {
  static_assert(sizeof(LidLvlPair) == sizeof(LidLvlPair::Id),
                "LidLvlPair Id sizeof assertion");
  bool ok = true;
  for (const Int lid : {0, 11, 128, (1 << 22) - 1})
    for (const Int lvl : {0, 1, 43, 256, 1023}) {
      LidLvlPair p(lid, lvl);
      if (p.lid() != lid || p.lvl() != lvl) ok = false;
    }
  return ok;
}

static void assert_sizes (int nentity, int nlvl) {
  assert(nlvl <= (1 << LidLvlPair::sz) - 1);
  assert(nentity <= (1 << (32 - LidLvlPair::sz)) - 1);
}

struct ExchangeMaskedData {
  const Distribution::Data& dd;
  const Array2D<const char>& mask;
  const Array3D<Real>* const fs;
  LidLvlPair* idxs;
  uint* ptrs;
  uint nidxs;
  
  ExchangeMaskedData (const Distribution::Data& dd_, const Array2D<const char>& mask_,
                      const Array3D<Real>* const fs_)
    : dd(dd_), mask(mask_), fs(fs_)
  {}
};

struct Comm::Impl {
  typedef uint32_t uint;

  Distribution::ConstPtr d;
  MpasMesh::ConstPtr m;
  std::vector<Real> srbuf, rrbuf;
  std::vector<LidLvlPair::Id> ribuf;
  std::vector<uint> em_iwrk;
  std::vector<mpi::Request> rreq, sreq, rreq1, sreq1;

  // Metadata for exchange_and_accumulate.
  struct EaaMetaData {
    Int nlvl;
    std::vector<size_t> rbuf_os;
    EaaMetaData () : nlvl(0) {}
  };
  EaaMetaData eaamd;

  Impl (Distribution::ConstPtr d_, MpasMesh::ConstPtr m_) : d(d_), m(m_) {}

  void exchange_full(const mpi::Parallel& p, const Aci1& pids,
                     const Distribution::Data& dd, const int nf, const int* ns,
                     const int nlvl, const Array3D<Real>* const fs);
  void exchange_masked(const mpi::Parallel& p, const Aci1& pids, const int nemd,
                       ExchangeMaskedData* const emd, const char maskval,
                       const int nf, const int* ns, const int nlvl);
  template <typename AccumT>
  void exchange_and_accumulate(const mpi::Parallel& p, const Aci1& pids,
                               const Distribution::Data& dd, const int nlvl,
                               const Array2D<AccumT>& f);

  bool test(const Distribution::Data& dd);
  bool test_exchange_and_accumulate(const Distribution::Data& dd);
  bool test_exchange_and_accumulate(const Distribution::Data& dd, const int fill);
  bool test();
};

void Comm::Impl
::exchange_full (const mpi::Parallel& p, const Aci1& pids,
                 const Distribution::Data& dd, const int nf, const int* ns,
                 const int nlvl, const Array3D<Real>* const fs) {
  if (p.size() == 1) return;

  const Int tag = 79;
  const Int ncown = dd.nown, npid = pids.m();
  const auto& pidptrs = dd.pidptrs;
  const auto& olids_ptrs = dd.olids_ptrs;
  const auto& olids = dd.olids;

  assert(pidptrs(0) == 0); assert(pidptrs(1) == ncown);
  assert(olids_ptrs(0) == 0); assert(olids_ptrs(1) == 0);
  assert(pidptrs(npid) == ncown + dd.rlids.m());
  assert(olids_ptrs(npid) == olids.m());

  Int ntot_subfield = 0;
  for (Int i = 0; i < nf; ++i) {
    ntot_subfield += ns[i];
    assert(fs[i].m() >= ns[i]);
    assert(fs[i].n() >= nlvl);
  }
  const Int nscal_per_col = nlvl*ntot_subfield;

  srbuf.resize(nscal_per_col*olids.m());
  rrbuf.resize(nscal_per_col*dd.rlids.m());
  rreq.resize(npid);
  sreq.resize(npid);

  // Set up recvs.
  for (Int pi = 1; pi < npid; ++pi) {
    const Int pid = pids(pi);
    assert(size_t(nscal_per_col)*(pidptrs(pi) - ncown) +
           nscal_per_col*(pidptrs(pi+1) - pidptrs(pi)) <= rrbuf.size());
    mpi::irecv(p, &rrbuf[nscal_per_col*(pidptrs(pi) - ncown)],
               nscal_per_col*(pidptrs(pi+1) - pidptrs(pi)),
               pid, tag, &rreq[pi]);
  }

  // Write send data to buf and send.
  for (Int pi = 1; pi < npid; ++pi) {
    const Int pid = pids(pi);
    const Int nent = olids_ptrs(pi+1) - olids_ptrs(pi);
    const size_t os = olids_ptrs(pi);
    Int fi_os = nscal_per_col*os;
    for (Int fi = 0; fi < nf; ++fi) {
      const Int nsubfield = ns[fi];
      const auto& f = fs[fi];
      for (Int i = 0; i < nent; ++i) {
        const Int enti = olids(os + i);
        for (Int li = 0; li < nlvl; ++li)
          for (Int qi = 0; qi < nsubfield; ++qi)
            srbuf[fi_os + (i*nlvl + li)*nsubfield + qi] = f(qi,li,enti);
      }
      fi_os += nent*nlvl*nsubfield;
    }
    assert(nscal_per_col*os + nscal_per_col*nent <= srbuf.size());
    mpi::isend(p, &srbuf[nscal_per_col*os], nscal_per_col*nent,
               pid, tag, &sreq[pi]);
  }

  // Wait on recvs and write to fields.
  for (Int i = 0; i < npid-1; ++i) {
    Int pi;
    mpi::waitany(npid-1, &rreq[1], &pi);
    ++pi;
    const Int os = pidptrs(pi), nent = pidptrs(pi+1) - pidptrs(pi);
    Int fi_os = nscal_per_col*(os - ncown);
    for (Int fi = 0; fi < nf; ++fi) {
      const Int nsubfield = ns[fi];
      auto& f = fs[fi];
      for (Int i = 0; i < nent; ++i) {
        for (Int li = 0; li < nlvl; ++li)
          for (Int qi = 0; qi < nsubfield; ++qi)
            f(qi,li,os+i) = rrbuf[fi_os + (i*nlvl + li)*nsubfield + qi];
      }
      fi_os += nent*nlvl*nsubfield;
    }
  }

  // Wait on sends.
  mpi::waitall(npid-1, &sreq[1]);
}

// Sorted list of src halo (LID,level) pairs to request. We actually record
// LID-nown.
static uint get_sorted_list_of_pairs_to_request (
  const Distribution::Data& dd, const int npid, const int nlvl,
  const Array2D<const char>& mask, const char maskval,
  // Pointers into idxs segmenting idxs into remotes.
  LidLvlPair* const idxs,
  // Number of halo (LID,level) pairs for which I need data.
  uint* const ptrs)
{
  const auto nown = dd.nown, nhalo = dd.nhalo;
  uint nidxs = 0;
  for (Int i = 0; i < nhalo; ++i)
    for (Int k = 0; k < nlvl; ++k)
      if (mask(k,i) & maskval)
        idxs[nidxs++].set(i, k);
  ptrs[0] = 0; // follow the convention that the first slot is for myself and unused
  ptrs[1] = 0;
  Int j = 2;
  for (uint i = 0; i < nidxs; ) {
    if (nown + idxs[i].lid() < dd.pidptrs(j)) {
      ++i;
    } else {
      ptrs[j] = i;
      ++j;
      assert(j <= npid);
      ptrs[j] = ptrs[j-1];
    }
  }
  for ( ; j <= npid; ++j)
    ptrs[j] = nidxs;
#ifndef NDEBUG
  for (Int i = 0; i < npid; ++i) assert(ptrs[i] <= ptrs[i+1]);
  assert(ptrs[npid] == nidxs);
#endif
  return nidxs;
}

// Called "upwind MPI" in some contexts.
void Comm::Impl
::exchange_masked (const mpi::Parallel& p, const Aci1& pids, const int nemd,
                   ExchangeMaskedData* const emd, const char maskval,
                   const int nf, const int* ns, const int nlvl) {
  if (p.size() == 1) return;

  const Int tag = 29;
  const auto npid = pids.m();

  size_t n = 2*npid;
  for (int ie = 0; ie < nemd; ++ie) {
    const auto& dd = emd[ie].dd;
    assert_sizes(dd.nown + dd.nhalo, nlvl);
    n += size_t(dd.nhalo)*nlvl + npid + 1;
  }
  if (em_iwrk.size() < n) em_iwrk.resize(n);

  if ((Int) rreq .size() < npid) { rreq .resize(npid); sreq .resize(npid); }
  if ((Int) rreq1.size() < npid) { rreq1.resize(npid); sreq1.resize(npid); }

  size_t niwrk = 0;
  auto* const tot_ptrs = em_iwrk.data() + niwrk;
  niwrk += npid + 1;
  for (int pi = 0; pi <= npid; ++pi) tot_ptrs[pi] = 0;
  for (int ie = 0; ie < nemd; ++ie) {
    auto& e = emd[ie];
    const auto& dd = e.dd;
    e.idxs = reinterpret_cast<LidLvlPair*>(em_iwrk.data() + niwrk);
    niwrk += dd.nhalo*nlvl;
    e.ptrs = em_iwrk.data() + niwrk;
    niwrk += npid + 1;
    e.nidxs = get_sorted_list_of_pairs_to_request(dd, npid, nlvl, e.mask,
                                                  maskval, e.idxs, e.ptrs);
    for (int pi = 0; pi <= npid; ++pi) tot_ptrs[pi] += e.ptrs[pi];
  }

  Int ntot_subfield = 0;
  for (Int i = 0; i < nf; ++i) {
    ntot_subfield += ns[i];
    for (int ie = 0; ie < nemd; ++ie) {
      auto& e = emd[ie];
      assert(e.fs[i].m() >= ns[i]);
      assert(e.fs[i].n() >= nlvl );
      assert(e.fs[i].o() >= e.dd.nhalo);
    }
  }

  // Set up recvs of remotes' requests for my LIDs. A request starts with a
  // count of (LID,lvl) pairs, then lists the pairs.
  n = 0;
  for (int ie = 0; ie < nemd; ++ie)
    n += size_t(nlvl)*emd[ie].dd.olids_ptrs(npid) + npid - 1;
  /*if (ribuf.size() < n)*/ ribuf.resize(n);
  for (Int pi = 1; pi < npid; ++pi) {
    const Int pid = pids(pi);
    size_t os = nemd*(pi - 1), sz = nemd; // os and size for max possible count
    for (int ie = 0; ie < nemd; ++ie) {
      const auto& dd = emd[ie].dd;
      os += nlvl*dd.olids_ptrs(pi);
      sz += nlvl*(dd.olids_ptrs(pi+1) - dd.olids_ptrs(pi));
    }
    assert(os + sz <= ribuf.size());
    mpi::irecv(p, &ribuf[os], sz, pid, tag, &rreq[pi]);
  }
    
  // Send my requests to remotes, translating my LIDs to their LIDs using rlids.
  {
    n = 0;
    for (int ie = 0; ie < nemd; ++ie)
      n += emd[ie].nidxs + npid - 1;
    n = (n*sizeof(LidLvlPair::Id) + sizeof(Real) - 1)/sizeof(Real);
    if (srbuf.size() < n) srbuf.resize(n);
    auto* const sibuf = reinterpret_cast<LidLvlPair::Id*>(srbuf.data());
    size_t cnt = 0;
    for (Int pi = 1; pi < npid; ++pi) {
      const Int pid = pids(pi);
      auto* const sbuf = &sibuf[tot_ptrs[pi] + nemd*(pi - 1)];
      size_t os = 0;
      for (int ie = 0; ie < nemd; ++ie) {
        const auto& e = emd[ie];
        const auto& dd = e.dd;
        const auto ptrs = e.ptrs;
        const Int npair = ptrs[pi+1] - ptrs[pi];
        sbuf[os++] = npair;
        const auto* const pairs = &e.idxs[ptrs[pi]];
        auto* const rpairs = reinterpret_cast<LidLvlPair*>(&sbuf[os]);
        os += npair;
        const auto& lidmap = dd.rlids; // my LID - nown -> remote LID
        ompparfor for (Int j = 0; j < npair; ++j)
          rpairs[j].set(lidmap[pairs[j].lid()], pairs[j].lvl());
      }
      cnt += os;
      assert(os == tot_ptrs[pi+1] - tot_ptrs[pi] + nemd);
      assert(cnt*sizeof(LidLvlPair::Id) <= n*sizeof(Real));
      assert(size_t(tot_ptrs[pi]) + nemd*(pi - 1) + os <=
             srbuf.size()*((sizeof(Real) + sizeof(LidLvlPair::Id) - 1)/
                           sizeof(LidLvlPair::Id)));
      mpi::isend(p, sbuf, os, pid, tag, &sreq[pi]);
    }
  }

  // Set up recvs of remotes' responses to my requests.
  n = 0;
  for (int ie = 0; ie < nemd; ++ie) n += emd[ie].nidxs;
  n *= ntot_subfield;
  if (rrbuf.size() < n) rrbuf.resize(n);
  auto* const nz_remotes = em_iwrk.data() + niwrk;
  niwrk += npid - 1;
  assert(niwrk <= em_iwrk.size());
  assert(em_iwrk.data() + niwrk <= nz_remotes + npid - 1);
  Int nnz_remotes = 0;
  for (Int pi = 1; pi < npid; ++pi) {
    const Int npair = tot_ptrs[pi+1] - tot_ptrs[pi];
    if (npair == 0) continue; // skip message if no data
    const Int pid = pids(pi);
    assert(size_t(ntot_subfield)*tot_ptrs[pi] + ntot_subfield*npair <=
           rrbuf.size());
    mpi::irecv(p, &rrbuf[ntot_subfield*tot_ptrs[pi]], ntot_subfield*npair,
               pid, tag+1, &rreq1[nnz_remotes]);
    nz_remotes[nnz_remotes++] = pi;
  }
  assert(nnz_remotes < npid);

  // Wait on remotes' requests.
  mpi::waitall(npid-1, &rreq[1]);

  // Compute send buf size. We did waitall rather than waitany so that all
  // counts would be available for this calculation.
  {
    Int nsend = 0;
    for (Int pi = 1; pi < npid; ++pi) {
      size_t os = nemd*(pi - 1);
      for (int ie = 0; ie < nemd; ++ie)
        os += nlvl*emd[ie].dd.olids_ptrs(pi);
      for (int ie = 0; ie < nemd; ++ie) {
        const auto& dd = emd[ie].dd;
        const Int npair = ribuf[os];
        os += npair + 1;
        assert(npair >= 0);
        assert(npair <= dd.nown*nlvl);
        nsend += npair;
      }
    }
    // Clean up the send requests so the send buffer is available.
    mpi::waitall(npid-1, &sreq[1]);
    n = size_t(ntot_subfield)*nsend;
    if (srbuf.size() < n) srbuf.resize(n);
  }

  // Fullfill remotes' requests.
  Int sn1 = 0;
  {
    size_t ros = 0;
    for (Int pi = 1; pi < npid; ++pi) {
      size_t ios = nemd*(pi - 1);
      for (int ie = 0; ie < nemd; ++ie)
        ios += nlvl*emd[ie].dd.olids_ptrs(pi);
      Real* const rbuf = &srbuf[ros];
      size_t cnt = 0;
      for (int ie = 0; ie < nemd; ++ie) {
        const auto& e = emd[ie];
        const uint npair = ribuf[ios];
        assert(npair <= uint(nlvl*e.dd.nown));
        ios++;
        if (npair == 0) continue;
        const auto* const pairs = reinterpret_cast<const LidLvlPair*>(&ribuf[ios]);
        ios += npair;
        Real* const rbuf_ie = rbuf + cnt;
        for (Int fi = 0, fi_os = 0; fi < nf; ++fi) {
          const auto nsubfield = ns[fi];
          const auto& f = e.fs[fi];
          ompparfor for (uint i = 0; i < npair; ++i) {
            const Int enti = pairs[i].lid(), li = pairs[i].lvl();
            assert(enti >= 0 && enti < e.dd.nown);
            for (Int qi = 0; qi < nsubfield; ++qi)
              rbuf_ie[fi_os + i*nsubfield + qi] = f(qi,li,enti);
          }
          fi_os += npair*nsubfield;
        }
        cnt += npair*ntot_subfield;
      }
      if (cnt == 0) continue;
      mpi::isend(p, rbuf, cnt, pids(pi), tag+1, &sreq1[sn1++]);
      ros += cnt;
    }
  }

  // Receive the data I requested.
  for (Int i = 0; i < nnz_remotes; ++i) {
    Int idx;
    mpi::waitany(nnz_remotes, rreq1.data(), &idx);
    const auto pi = nz_remotes[idx];
    const auto* const rbuf = &rrbuf[ntot_subfield*tot_ptrs[pi]];
    uint fi_os = 0;
    for (int ie = 0; ie < nemd; ++ie) {
      const auto& e = emd[ie];
      const auto* const pairs = &e.idxs[e.ptrs[pi]];
      const auto npair = e.ptrs[pi+1] - e.ptrs[pi];
      assert(npair >= 0);
      if (npair == 0) continue;
      for (Int fi = 0; fi < nf; ++fi) {
        const auto nsubfield = ns[fi];
        auto& f = e.fs[fi];
        ompparfor for (uint i = 0; i < npair; ++i) {
          const auto enti = e.dd.nown + pairs[i].lid(), li = pairs[i].lvl();
          for (Int qi = 0; qi < nsubfield; ++qi)
            f(qi,li,enti) = rbuf[fi_os + i*nsubfield + qi];
        }
        fi_os += npair*nsubfield;
      }
    }
  }

  // Clean up the send requests.
  mpi::waitall(sn1, sreq1.data());
}

static Int eaa_calc_nreal_metadata (const Int n) {
  // 1 in n+1 is for the integer at the start, given that sizeof(Int) <=
  // sizeof(LidLvlPair).
  static_assert(sizeof(Int) <= sizeof(LidLvlPair),
                "LidLvlPair Int sizeof assertion");
  return ((n+1)*sizeof(LidLvlPair) + sizeof(Real) - 1)/sizeof(Real);
}

typedef OrderInvariantAccumTraits::LocalAtomic LocalAtomic;

template <typename T> struct RepTraits;
template <> struct RepTraits<Real> {
  typedef Real type;
  enum : int { nreal = 1 };
};
template <> struct RepTraits<LocalAtomic> {
  typedef LocalAtomic::IRep type;
  enum : int { nreal = OrderInvariantAccumTraits::local_atomic_size_nreal };
};

inline bool is_zero (const Real& a) { return a == 0; }
inline Real get_rep (const Real& a) { return a; }
inline LocalAtomic::IRep get_rep (const LocalAtomic& a) { return a.get_irep(); }

template <typename T> inline void
accum_slot (const Array2D<T>& f, const int lid, const int lvl,
            const typename T::IRep& term) {
  auto& f_lvl_lid = f(lvl,lid);
  f_lvl_lid += T(f_lvl_lid.get_maxmag(), term);
}

static inline void
accum_slot (const Ar2& f, const int lid, const int lvl, const Real term) {
  auto& f_lvl_lid = f(lvl,lid);
  ompatomic f_lvl_lid += term;
}

// Works only for cell values because edge ownership is more complicated. This
// can be considered "downwind MPI" and requires just one comm round instead of
// two, but with the same comm volume. The (lid,lvl) list is in the same message
// as the bulk data.
template <typename AccumT>
void Comm::Impl
::exchange_and_accumulate (const mpi::Parallel& p, const Aci1& pids,
                           const Distribution::Data& dd, const int nlvl,
                           const Array2D<AccumT>& f) {
  if (p.size() == 1) return;

  typedef typename RepTraits<AccumT>::type Rep;
  static constexpr int rep_nreal = RepTraits<AccumT>::nreal;

  // Only one comm round is needed. I send pids(i) the following message
  // covering entities pidptrs(i):pidptrs(i+1)-1:
  //   (n LidLvlPair*n Real*n),
  // with padding as needed for alignment.
  const Int tag = 17;
  const auto nown = dd.nown, nhalo = dd.nhalo, npid = pids.m();
  assert_sizes(nown + nhalo, nlvl);
  const auto& pidptrs = dd.pidptrs;
  const auto& olids_ptrs = dd.olids_ptrs;
  const auto& rlids = dd.rlids;
  const Int nscal_per_col = nlvl;
  size_t bos;

  if (eaamd.nlvl != nlvl) {
    eaamd.rbuf_os.resize(npid+1, 0);
    bos = 0;
    for (Int pi = 1; pi < npid; ++pi) {
      eaamd.rbuf_os[pi] = bos;
      const Int nbulk = nscal_per_col*(olids_ptrs(pi+1) - olids_ptrs(pi));
      bos += eaa_calc_nreal_metadata(nbulk) + nbulk*rep_nreal;
    }
    eaamd.rbuf_os[npid] = bos;
    if (rrbuf.size() < bos) rrbuf.resize(bos);
    bos = 0;
    for (Int pi = 1; pi < npid; ++pi) {
      const Int nbulk = nscal_per_col*(pidptrs(pi+1) - pidptrs(pi));
      bos += eaa_calc_nreal_metadata(nbulk) + nbulk*rep_nreal;
    }
    if (srbuf.size() < bos) srbuf.resize(bos);
    rreq.resize(npid);
    sreq.resize(npid);
    eaamd.nlvl = nlvl;
  }

  // Set up receives.
  for (Int pi = 1; pi < npid; ++pi)
    mpi::irecv(p, &rrbuf[eaamd.rbuf_os[pi]],
               eaamd.rbuf_os[pi+1] - eaamd.rbuf_os[pi],
               pids(pi), tag, &rreq[pi]);

  // Pack and send.
  bos = 0;
  for (Int pi = 1; pi < npid; ++pi) {
    Real* const rptr = srbuf.data() + bos;
    const Int nent = pidptrs(pi+1) - pidptrs(pi);
    const Int nbulk = nscal_per_col*nent;
    bos += eaa_calc_nreal_metadata(nbulk) + nbulk*rep_nreal;
    // Count non-0s entries and record (rlid,lvl).
    void* const ptr = rptr;
    const Int os = pidptrs(pi);
    auto* const siptr = reinterpret_cast<LidLvlPair*>(ptr) + 1;
    Int n = 0;
    for (Int i = 0; i < nent; ++i)
      for (Int li = 0; li < nlvl; ++li) {
        if (is_zero(f(li,os+i))) continue;
        siptr[n++].set(rlids(os-nown+i), li);
        assert(siptr[n-1].lid() == rlids(os-nown+i));
        assert(siptr[n-1].lvl() == li);
      }
    *reinterpret_cast<Int*>(ptr) = n;
    // Record bulk data.
    const Int ros = eaa_calc_nreal_metadata(n);
    auto* const srptr = reinterpret_cast<Rep*>(rptr + ros);
    n = 0;
    for (Int i = 0; i < nent; ++i)
      for (Int li = 0; li < nlvl; ++li) {
        const auto& r = f(li,os+i);
        if (is_zero(r)) continue;
        srptr[n++] = get_rep(r);
      }
    mpi::isend(p, rptr, ros + n*rep_nreal, pids(pi), tag, &sreq[pi]);
  }

  // Accumulate receives.
  for (Int pi = 1; pi < npid; ++pi) {
    Int idx;
    mpi::waitany(npid-1, &rreq[1], &idx);
    bos = eaamd.rbuf_os[idx+1];
    Real* const rptr = rrbuf.data() + bos;
    void* const ptr = rptr;
    const Int n = *reinterpret_cast<Int*>(ptr);
    const auto* const riptr = reinterpret_cast<LidLvlPair*>(ptr) + 1;
    const Int ros = eaa_calc_nreal_metadata(n);
    const auto* const rrptr = reinterpret_cast<Rep*>(rptr + ros);
    ompparfor for (Int i = 0; i < n; ++i) {
      const auto& p = riptr[i];
      assert(p.lid() >= 0); assert(p.lid() < nown);
      assert(p.lvl() >= 0); assert(p.lvl() < nlvl);
      accum_slot(f, p.lid(), p.lvl(), rrptr[i]);
    }
  }

  mpi::waitall(npid-1, &sreq[1]);
}

bool Comm::Impl
::test_exchange_and_accumulate (const Distribution::Data& dd, const int fill) {
  const Real eps = std::numeric_limits<Real>::epsilon();
  bool ok = true;
  const Int npid = d->pids.m(), nlvl = 3;
  const auto& p = *d->p;
  const Array2D<Real> f(nlvl, dd.gids.m()), f_t(nlvl, dd.nown);
  f.set(0); f_t.set(0);
  for (Int i = 0; i < dd.gids.m(); ++i)
    for (Int li = 0; li < nlvl; ++li) {
      if (fill == 1) {
        if ((dd.gids(i) + li) % (p.rank() + 4) == 0)
          f(li,i) = (11*dd.gids(i) + 7*li) % 7919;
      } else {
        f(li,i) = square(p.rank()) + li;
      }
      if (i < dd.nown) f_t(li,i) = f(li,i);
    }
  exchange_and_accumulate(p, d->pids, dd, nlvl, f);
  for (Int pid = 0; pid < npid; ++pid) {
    const Int rank = d->pids(pid);
    for (Int j = dd.olids_ptrs(pid); j < dd.olids_ptrs(pid+1); ++j) {
      const Int i = dd.olids(j);
      for (Int li = 0; li < nlvl; ++li) {
        if (fill == 1) {
          if ((dd.gids(i) + li) % (rank + 4) == 0)
            f_t(li,i) += (11*dd.gids(i) + 7*li) % 7919;
        } else {
          f_t(li,i) += square(rank) + li;
        }
      }
    }
  }
  for (Int i = 0; i < dd.nown; ++i)
    for (Int li = 0; li < nlvl; ++li)
      // Order of ops differences in forming f and f_t can result in some error.
      //todo-pebfb Once I have invariant summation, use it here.
      if (reldif(f(li,i), f_t(li,i)) > 10*eps)
        ok = false;
  return ok;
}

bool Comm::Impl::test_exchange_and_accumulate (const Distribution::Data& dd) {
  const bool t1 = test_exchange_and_accumulate(dd, 1);
  const bool t2 = test_exchange_and_accumulate(dd, 2);
  return t1 && t2;
}

bool Comm::Impl::test (const Distribution::Data& dd) {
  static const Int nf = 2, nlvl = 3;
  const Int ns[nf] = {1, 4};
  Array3D<Real> fs[nf];
  const auto fill = [&] () {
    for (Int j = 0; j < nf; ++j) {
      auto& f = fs[j];
      f.reset(ns[j], nlvl, dd.gids.m());
      f.set(1e30);
      for (Int i = 0; i < dd.nown; ++i)
        for (Int li = 0; li < nlvl; ++li)
          for (Int qi = 0; qi < ns[j]; ++qi)
            f(qi,li,i) = dd.gids(i) + qi + 0.1*li;
    }
  };
  bool ok = LidLvlPair::test();
  const auto test = [&] (const Distribution::Data& dd, const int os,
                         const int incr[nlvl]) {
    int tested = 0;
    for (Int j = 0; j < nf; ++j) {
      const auto& f = fs[j];
      for (Int li = 0; li < nlvl; ++li)
        for (Int qi = 0; qi < ns[j]; ++qi)
          for (Int i = os; i < dd.gids.m(); i += incr[li]) {
            ++tested;
            if (f(qi,li,i) != Real(dd.gids(i) + qi + 0.1*li))
              ok = false;
          }
    }
    return tested;
  };
  {
    static const int incr[nlvl] = {1, 1, 1};
    fill();
    exchange_full(*d->p, d->pids, dd, nf, ns, nlvl, fs);
    test(dd, 0, incr);
  }
  {
    static const int incr[nlvl] = {1, 1, 1};
    fill();
    const Array2D<char> mask(nlvl, dd.nhalo);
    mask.set(1);
    ExchangeMaskedData emd(dd, mask, fs);
    exchange_masked(*d->p, d->pids, 1, &emd, 1, nf, ns, nlvl);
    test(dd, 0, incr);
  }
  {
    static const int incr[nlvl] = {1, 2, 7};
    fill();
    const Array2D<char> mask(nlvl, dd.nhalo);
    mask.set(0);
    for (Int li = 0; li < nlvl; ++li)
      for (Int i = 0; i < dd.nhalo; i += incr[li])
        mask(li,i) = 1;
    ExchangeMaskedData emd(dd, mask, fs);
    exchange_masked(*d->p, d->pids, 1, &emd, 1, nf, ns, nlvl);
    test(dd, dd.nown, incr);
  }
  {
    static const int incr[nlvl] = {1, 2, 7};
    fill();
    const Array2D<char> mask(nlvl, dd.nhalo);
    mask.set(0);
    for (Int li = 0; li < nlvl; ++li)
      for (Int i = 0; i < dd.nhalo; i += incr[li])
        mask(li,i) = 1;
    ExchangeMaskedData
      emd[] = {ExchangeMaskedData(dd, mask, fs),
               ExchangeMaskedData(dd, mask, fs)};
    exchange_masked(*d->p, d->pids, 2, emd, 1, nf, ns, nlvl);
    test(dd, dd.nown, incr);
  }
  return ok;
}

bool Comm::Impl::test () {
  const bool c = test(d->c), e = test(d->e);
  const bool eaa = test_exchange_and_accumulate(d->c);
  return c && e && eaa;
}

Comm::Comm (Distribution::ConstPtr d, MpasMesh::ConstPtr m) {
  p = std::make_shared<Impl>(d, m);
  p->d = d;
}

int Comm::rank () const { return p->d->p->rank(); }
int Comm::nrank () const { return p->d->p->size(); }
bool Comm::amroot () const { return p->d->p->amroot(); }
int Comm::szhalo () const { return p->d->szhalo; }

const Array1D<const char>& Comm::cdhalo () const {
  static const Array1D<const char> cdhalo = p->d->cdhalo;
  return cdhalo;
}

Distribution::ConstPtr Comm::get_distribution () const { return p->d; }

void Comm
::exchange_cell_full (const int nf, const int* ns, const int nlvl,
                      const Array3D<Real>* const fs) {
  p->exchange_full(*p->d->p, p->d->pids, p->d->c, nf, ns, nlvl, fs);
}

void Comm
::exchange_cell_masked (const Array2D<const char>& mask, const char maskval,
                        const int nf, const int* ns, const int nlvl,
                        const Array3D<Real>* const fs) {
  ExchangeMaskedData emd(p->d->c, mask, fs);
  p->exchange_masked(*p->d->p, p->d->pids, 1, &emd, maskval, nf, ns, nlvl);
}

void Comm::exchange_and_accumulate_cell (const int nlvl, const Array2D<Real>& f) {
  p->exchange_and_accumulate(*p->d->p, p->d->pids, p->d->c, nlvl, f);
}

void Comm::exchange_and_accumulate_cell (
  const int nlvl, const Array2D<LocalAtomic>& f)
{ p->exchange_and_accumulate(*p->d->p, p->d->pids, p->d->c, nlvl, f); }

void Comm
::exchange_edge_full (const int nf, const int* ns, const int nlvl,
                      const Array3D<Real>* const fs) {
  p->exchange_full(*p->d->p, p->d->pids, p->d->e, nf, ns, nlvl, fs);
}

void Comm
::exchange_edge_masked  (const Array2D<const char>& mask, const char maskval,
                         const int nf, const int* ns, const int nlvl,
                         const Array3D<Real>* const fs) {
  ExchangeMaskedData emd(p->d->e, mask, fs);
  p->exchange_masked(*p->d->p, p->d->pids, 1, &emd, maskval, nf, ns, nlvl);
}

void Comm
::exchange_edge_and_cell_masked (
  const Array2D<const char>& edge_mask, const Array2D<const char>& cell_mask,
  const char maskval, const int nf, const int* ns, const int nlvl,
  const Array3D<Real>* const edge_fs, const Array3D<Real>* const cell_fs)
{
  ExchangeMaskedData
    emd[2] = {ExchangeMaskedData(p->d->e, edge_mask, edge_fs),
              ExchangeMaskedData(p->d->c, cell_mask, cell_fs)};
  p->exchange_masked(*p->d->p, p->d->pids, 2, emd, maskval, nf, ns, nlvl);
}

static MPI_Op get_mpi_op (const Comm::Op op)
{ return op == Comm::SUM ? MPI_SUM : op == Comm::MAX ? MPI_MAX : MPI_MIN; }

void Comm::all_reduce (const Op op, const Int n, const Real* l, Real* g)
{ mpi::all_reduce(*p->d->p, l, g, n, get_mpi_op(op)); }
void Comm::all_reduce (const Op op, const Int n, const Int* l, Int* g)
{ mpi::all_reduce(*p->d->p, l, g, n, get_mpi_op(op)); }
void Comm::all_reduce (const Op op, const Int n, const int64_t* l, int64_t* g)
{ mpi::all_reduce(*p->d->p, l, g, n, get_mpi_op(op)); }

static void
kbn_sum (void* invec, void* inoutvec, int* len, MPI_Datatype* datatype) {
  assert(*len % 2 == 0);
  const int n = *len/2;
  const auto* s = reinterpret_cast<const KahanBabushkaNeumaierSummer*>(invec);
  auto* d = reinterpret_cast<KahanBabushkaNeumaierSummer*>(inoutvec);
  for (int i = 0; i < n; ++i)
    d[i] += s[i];
}

void Comm
::all_reduce (const Op op, const Int n, const KahanBabushkaNeumaierSummer* l,
              KahanBabushkaNeumaierSummer* g) {
  assert(op == Comm::SUM); // for now; later, could support min, max
  typedef KahanBabushkaNeumaierSummer::scalar_type ST;
  const ST* s = reinterpret_cast<const ST*>(l);
  ST* r = reinterpret_cast<ST*>(g);
  mpi::Op mop(kbn_sum, true);
  mpi::all_reduce(*p->d->p, s, r, 2*n, mop);
}

static void
ld_sum (void* invec, void* inoutvec, int* len, MPI_Datatype* datatype) {
  assert(*len % LongDoubleSummer::nint == 0);
  const int n = *len/LongDoubleSummer::nint;
  const auto* s = reinterpret_cast<const LongDoubleSummer*>(invec);
  auto* d = reinterpret_cast<LongDoubleSummer*>(inoutvec);
  for (int i = 0; i < n; ++i)
    d[i] += s[i];
}

void Comm
::all_reduce (const Op op, const Int n, const LongDoubleSummer* l,
              LongDoubleSummer* g) {
  assert(op == Comm::SUM); // for now; later, could support min, max
  const int* s = reinterpret_cast<const int*>(l);
  int* r = reinterpret_cast<int*>(g);
  mpi::Op mop(ld_sum, true);
  mpi::all_reduce(*p->d->p, s, r, LongDoubleSummer::nint*n, mop);
}

template <typename Int>
static void
Int2_sum (void* invec, void* inoutvec, int* len, MPI_Datatype* datatype) {
  static constexpr int nint = (2*sizeof(Int))/sizeof(int);
  assert(*len % nint == 0);
  const int n = *len/nint;
  // double doesn't matter, as it's unused.
  typedef order_invariant_accum::RealInt2Rep<double,Int> IntRep;
  const auto* s = reinterpret_cast<const IntRep*>(invec);
  auto* d = reinterpret_cast<IntRep*>(inoutvec);
  for (int i = 0; i < n; ++i)
    d[i] += s[i];
}

void Comm
::all_reduce (const Op op, const Int n,
              const order_invariant_accum::IntRep<double,true,false>* l,
              order_invariant_accum::IntRep<double,true,false>* g) {
  assert(op == Comm::SUM); // for now; later, could support min, max
  const int* s = reinterpret_cast<const int*>(l);
  int* r = reinterpret_cast<int*>(g);
  typedef order_invariant_accum::IntRep<double,true,false>::Int IntT;
  static constexpr int nint = (2*sizeof(IntT))/sizeof(int);
  mpi::Op mop(Int2_sum<IntT>, true);
  mpi::all_reduce(*p->d->p, s, r, nint*n, mop);
}

static int rar_maxsum_nscalar;

// MPI_MAX on Real elements 0:rar_maxsum_nscalar-2, MPI_SUM on int element
// rar_maxsum_nscalar-1.
static void
rar_maxsum (void* invec, void* inoutvec, int* len, MPI_Datatype* datatype) {
  typedef typename RealTypeTraits<Real>::SameSizeInt RInt;
  static_assert(sizeof(RInt) == sizeof(Real), "reinterpret_cast");
  const int n = *len/rar_maxsum_nscalar, nscalar = rar_maxsum_nscalar - 1;
  const auto* s = reinterpret_cast<const Real*>(invec);
  auto* d = reinterpret_cast<Real*>(inoutvec);
  auto di = *reinterpret_cast<RInt*>(&d[nscalar]);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < nscalar; ++j)
      d[j] = std::max(d[j], s[i*rar_maxsum_nscalar + j]);
    di += *reinterpret_cast<const RInt*>(&s[i*rar_maxsum_nscalar + nscalar]);
  }
  const auto* const buf = reinterpret_cast<const Real*>(&di);
  d[nscalar] = *buf;
}

void Comm
::repro_all_reduce (const Op op, const int nterm, const Int nscalar,
                    const Ar1& wrk, const Real* l, Real* g) {
  assert(op == Comm::SUM);
  typedef order_invariant_accum::AccumType<Real,true,false>::type Accum;

  assert(wrk.m() >= nscalar*Accum::nreal + nscalar + 1);
  const auto os = nscalar*Accum::nreal;
  auto* const s = reinterpret_cast<Accum*>(wrk.data());
  
  {
    typedef typename RealTypeTraits<Real>::SameSizeInt RInt;
    Real* const lmaxmag = wrk.data() + os;
    Real* const gmaxmag = lmaxmag + nscalar + 1;
    calc_maxmag(nterm, nscalar, l, lmaxmag);
    RInt gnterm;
    rar_maxsum_nscalar = nscalar + 1;
    const RInt v = nterm;
    const auto* const buf = reinterpret_cast<const Real*>(&v);
    lmaxmag[nscalar] = *buf;
    mpi::Op mop(rar_maxsum, true);
    mpi::all_reduce(*p->d->p, lmaxmag, gmaxmag, rar_maxsum_nscalar, mop);
    gnterm = *reinterpret_cast<RInt*>(&gmaxmag[nscalar]);
    for (int j = 0; j < nscalar; ++j) {
      s[j].set_maxmag(gnterm*gmaxmag[j]);
      s[j].set_irep(0);
    }
  }

  Ar1 w1(wrk.m() - os, wrk.data() + os);
  Accum::sum(w1, l, nterm, s, nscalar, true);

  auto* const lirep = reinterpret_cast<Accum::IRep*>(wrk.data() + os);
  auto* const girep = lirep + nscalar;
  for (int j = 0; j < nscalar; ++j)
    lirep[j] = s[j].get_irep();

  all_reduce(op, nscalar, lirep, girep);

  for (int j = 0; j < nscalar; ++j)
    g[j] = Accum(s[j].get_maxmag(), girep[j])();
}

bool Comm::test () {
  // Run a bunch of tests in Impl.
  bool ok = p->test();
  { // Test KahanBabushkaNeumaierSummer reduction.
    static const int ntest = 5;
    const auto& par = *p->d->p;
    const int n = par.size(), rank = par.rank();
    const auto f = [&] (const int rank, const int t) {
      const int sign = t % 2 == 0 ? 1 : -1;
      return sign*std::sqrt(Real(t+1)/5 + rank);
    };
    KahanBabushkaNeumaierSummer lsum[ntest], gsum[ntest], tsum[ntest];
    for (int t = 0; t < ntest; ++t)
      lsum[t] += f(rank,t);
    all_reduce(Comm::SUM, ntest, lsum, gsum);
    for (int t = 0; t < ntest; ++t)
      for (int r = 0; r < n; ++r)
        tsum[t] += f(r,t);
    for (int t = 0; t < ntest; ++t)
      if (gsum[t]() != tsum[t]()) {
        ok = false;
        printf("%1d %23.16e %23.16e %1.3e\n",
               t, tsum[t](), gsum[t](), reldif(tsum[t](), gsum[t]()));
      }
  }
  return ok;
}

Comm::Ptr Comm::create (Distribution::ConstPtr d, MpasMesh::ConstPtr m) {
  return std::make_shared<Comm>(d, m);
}

const Int* get_cell_gids (const Comm& c) { return c.get_distribution()->c.gids.data(); }
const Int* get_edge_gids (const Comm& c) { return c.get_distribution()->e.gids.data(); }
