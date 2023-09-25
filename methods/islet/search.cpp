#ifndef SLMM_NP_GT_4
# define SLMM_NP_GT_4
#endif

#include "islet_tables.hpp"
#include "islet_npx.hpp"
#include "islet_maxeigcomp.hpp"
#include "islet_xnodes_metrics.hpp"
#include "islet_pum.hpp"
#include "islet_util.hpp"

#include <omp.h>
#include <array>
#include <algorithm>

class SearchAtom : public UserInterpMethod {
public:
  struct Input {
    enum Basis { gll, uniform, legendre, cheb };

    static const int np_max = 12;

    Int np, nmodregions;
    Basis basis;
    Int stabnp[np_max], staboffset[np_max];
    // Looking for at least this max |lambda| - 1.
    Real maxeigampm1;
    // Conclude the search when it succeeds with these parameters. The second is
    // the number of eigenvalues to sample in (0, 1/2].
    Int ne, neigdx;
    bool quiet, unittest;

    Input () { init(); }

  private:
    void init () {
      quiet = unittest = false;
      np = 6;
      basis = gll;
      nmodregions = 2;
      stabnp[0] = 6;
      stabnp[1] = 5;
      staboffset[0] = staboffset[1] = 0;
      maxeigampm1 = 1e-13;
      ne = 1000;
      neigdx = 1000;
    }
  };

  SearchAtom (const Input& iin) {
    reset(iin);
    if (in.unittest)
      std::cerr << (unittest() > 0 ? "FAIL" : "PASS")
                << ": SearchAtom::unittest.\n";
  }

  explicit SearchAtom () {}

  void reset (const Input& iin) {
    in = iin;
    switch (in.basis) {
    case Input::gll: x_gll = islet::get_x_gll_special(in.np); break;
    case Input::uniform: {
      static Real x[islet::np_max];
      x_gll = x;
      for (Int i = 0; i < in.np; ++i)
        x[i] = 2*(Real(i)/(in.np-1)) - 1;
    } break;
    case Input::legendre: x_gll = islet::get_x_gl(in.np); break;
    case Input::cheb: {
      static Real x[islet::np_max];
      x_gll = x;
      for (Int i = 0; i < in.np; ++i)
        x[i] = -std::cos(M_PI*Real(2*(i+1) - 1)/Real(2*in.np));
    } break;
    }
  }

  Real run () {
    return max_eig_amp.run(in.np, in.ne, in.neigdx, in.maxeigampm1,
                           in.quiet, this);
  }

  MaxEigComputer::Analysis calc_max_vals (const Int& nmu, const Int& ndx) {
    return max_eig_amp.calc_max_vals(nmu, ndx, in.np, this);
  }

  static int unittest () {
    int nerr = 0;
    Input in;
    SearchAtom sa(in);
    Real v0[32], v1[32];
    const Int np = 6;
    for (Int ix = 0, nx = 11; ix < nx; ++ix) {
      const auto x = -1 + (2.0*ix)/nx;
      npxstab<Real>::eval(np, x, v0);
      sa.eval(x, v1);
      for (Int j = 0; j < np; ++j)
        if (v0[j] != v1[j])
          ++nerr;
    }
    return nerr;
  }

  void eval (const Real& x, Real* const v) override {
    eval(in.np, in.nmodregions, x_gll,
         in.stabnp, in.staboffset,
         x, v);
  }

  Int get_np () const override { return in.np; }
  const Real* get_xnodes () const override { return x_gll; }

private:
  Input in;
  MaxEigComputer max_eig_amp;
  const Real* x_gll;
  std::vector<Real> x_gll_buf;

  static void eval (
    const Int& np, const Int& nreg, const Real* x_gll,
    const Int* const subnp, const Int* const os,
    const Real& x, Real* const v)
  {
    if (x > 0) {
      eval(np, nreg, x_gll, subnp, os, -x, v);
      for (int i = 0; i < np/2; ++i)
        std::swap(v[i], v[np-i-1]);
      return;
    }
    bool done = false;
    for (Int i = 0; i < nreg; ++i) {
      if (x > x_gll[i+1]) continue;
      assert(i == 0 || x >= x_gll[i]);
      assert( ! done);
      done = true;
      if (subnp[i] == np) {
        eval_lagrange_poly(x_gll, np, x, v);
      } else {
        std::fill(v, v + np, 0);
        eval_lagrange_poly(x_gll + os[i], subnp[i], x, v + os[i]);
      }
      break;
    }
    if ( ! done)
      eval_lagrange_poly(x_gll, np, x, v);
  }
};

static void calc_wts_metrics (const Int np, const Real* wt,
                              bool& all_pve_wts, Real& ratio) {
  all_pve_wts = true;
  for (Int i = 0; i < np; ++i) if (wt[i] <= 0) all_pve_wts = false;
  Real wtmin = 10, wtmax = -1;
  for (Int i = 0; i < np; ++i) wtmin = std::min(wtmin, wt[i]);
  for (Int i = 0; i < np; ++i) wtmax = std::max(wtmax, wt[i]);
  ratio = wtmax/wtmin;
}

static Real calc_pum_metric (UserInterpMethod& im, const bool threaded = true,
                             const Real stop_if_above = 1e3) {
  std::shared_ptr<UserInterpMethod> uim(&im, [] (UserInterpMethod*) {});
  const auto wrapper = std::make_shared<InterpMethod>(uim);
  pum::Options o;
  o.threaded = threaded;
  o.ntrial = 48;
  o.perturb = 0.01;
  Real pum_metric = 0;
  for (const Int mec_ne: {2, 4, 8, 16, 32, 64})
    for (const Int ne: {3, 5, 10}) {
      o.mec_ne = mec_ne;
      o.ne = ne;
      pum::PerturbedUniformMeshMetric pum(wrapper, o);
      const auto pum_metric_ne = pum.run(stop_if_above);
      // If we stopped, then we have no idea what the actual pum is, so return 1
      // to be safe.
      if (pum_metric > stop_if_above) return 1;
      pum_metric = std::max(pum_metric, pum_metric_ne);
    }
  return pum_metric;
}

// Restriction of nodal subset bases to an offset followed by adjacent nodes.
namespace find_offset_nodal_subset_bases {
static const int nps[] = {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16};

struct Basis {
  static const int M = islet::np_max;
  std::array<int,M> subnp, offst;
  int n;
};

bool same (const Basis& a, const Basis& b) {
  assert(a.n == b.n);
  bool s = true;
  for (int i = 0; i < a.n; ++i)
    if (a.subnp[i] != b.subnp[i] || a.offst[i] != b.offst[i]) {
      s = false;
      break;
    }
  return s;
}

static bool issymmetric (const int np, const int subnp, const int offst) {
  const int d = np - subnp;
  return ((d % 2 == 0) && offst == d/2);
}

static bool nodes_on_bdy (const SearchAtom::Input::Basis b) {
  return b == SearchAtom::Input::gll || b == SearchAtom::Input::uniform;
}

void recur (const int np, const int min_np, const int min_np_pos,
            std::vector<Basis>& good_bases, Basis b, const int pos,
            SearchAtom& sa, const SearchAtom::Input::Basis basis,
            MetricsTracker& mt) {
  if (pos == -1) {
    // Don't check if we've seen this basis is already good.
    bool already = false;
    for (const auto& gb : good_bases)
      if (same(gb, b)) {
        already = true;
        break;
      }
    if (already) return;

    SearchAtom::Input in;
    in.quiet = true;
    in.basis = basis;
    in.np = np;
    in.nmodregions = b.n;
    for (int i = 0; i < b.n; ++i) in.stabnp[i] = b.subnp[i];
    for (int i = 0; i < b.n; ++i) in.staboffset[i] = b.offst[i];
    in.ne = 500;
    in.neigdx = 500;
    sa.reset(in);

    Nodes nodes(np, nodes_on_bdy(basis)); {
      assert(nodes.get_nh() == b.n);
      std::vector<Int> ns;
      ns.reserve(np);
      for (Int ireg = 0; ireg < nodes.get_nh(); ++ireg) {
        ns.resize(in.stabnp[ireg]);
        for (Int i = 0; i < in.stabnp[ireg]; ++i)
          ns[i] = in.staboffset[ireg] + i;
        nodes.set(ireg, ns);
      }
    }
    Real xnodes_metric[3];
    calc_xnodes_metrics(nodes, sa.get_xnodes(), xnodes_metric);
    if ( ! mt.acceptable_metrics(nodes, sa.get_xnodes(), xnodes_metric)) return;
    const Real pum_to_accept = mt.pum_to_accept(nodes, sa.get_xnodes(),
                                                xnodes_metric);
    bool all_pve_wts = false;
    Real wtr = 0; {
      Real wt[islet::np_max];
      calc_weights(nodes, sa.get_xnodes(), wt);
      calc_wts_metrics(np, wt, all_pve_wts, wtr);
    }
    if ( ! all_pve_wts) return;

    // Run the potentially expensive analysis.
    Real maxeigampm1 = sa.run();
    Real maxcondV = 0, maxdefub = 0;
    Real pum_metric = 0;
    if (maxeigampm1 <= in.maxeigampm1) {
      pum_metric = calc_pum_metric(sa, true, pum_to_accept);
      if ( ! mt.would_update(xnodes_metric, pum_metric)) return;
      const auto max_vals = sa.calc_max_vals(1111, 1111);
      maxeigampm1 = max_vals.max_eig_amp - 1;
      maxcondV = max_vals.max_condv;
      maxdefub = max_vals.max_defect_ub;
    }
    if (maxeigampm1 <= in.maxeigampm1) {
      printf("meam1 %1.1e mcV %1.1e mdef %1.1e w>0 %d wtr %9.2e "
             "npm %9.2e %9.2e %9.2e pum %9.2e",
             maxeigampm1, maxcondV, maxdefub, all_pve_wts, wtr,
             xnodes_metric[0], xnodes_metric[1], xnodes_metric[2],
             pum_metric);
      printf(" | np %2d subnp", np);
      for (int i = 0; i < b.n; ++i) printf(" %d", b.subnp[i]);
      printf(" offst");
      for (int i = 0; i < b.n; ++i) printf(" %d", b.offst[i]);
      printf("\n");
      fflush(stdout);
      good_bases.push_back(b);
      mt.update(xnodes_metric, pum_metric);
    }
    return;
  }
  // Set up a basis. Avoid some, but not all, redundant trials by making one
  // slot have np = min_np.
  const bool middle = np % 2 == 0 && pos == np/2 - 1;
  for (
#if 0
    int subnp = pos == min_np_pos ? min_np : np;
    subnp >= min_np;
    --subnp
#else
    int subnp = min_np;
    subnp <= (pos == min_np_pos ? min_np : np);
    ++subnp
#endif
       )
    for (int offst = std::max(0, pos - subnp + 2);
         offst <= std::min(pos, np - subnp);
         ++offst) {
      if (middle && ! issymmetric(np, subnp, offst)) continue;
      b.subnp[pos] = subnp;
      b.offst[pos] = offst;
      recur(np, min_np, min_np_pos, good_bases, b, pos-1, sa, basis, mt);
    }
}

static Int run (const int np, MetricsTracker::Ptr mt = nullptr,
                const SearchAtom::Input::Basis basis = SearchAtom::Input::gll) {
  if ( ! mt) mt = std::make_shared<MetricsTracker>(np);
  const int min_ooa = 2; //np/2;
  const int min_np_lim = min_ooa;
  SearchAtom sa;
  std::vector<Basis> good_bases;
  printf("np %2d\n", np);
  Int good = 0, max_good_np = -1;
  // Search from highest- to lowest-OOA bases.
  for (int min_np = np; min_np >= min_np_lim; --min_np) {
    fprintf(stdout, "min_np %2d\n", min_np); fflush(stdout);
    Basis b;
    b.n = nodes_on_bdy(basis) ? np/2 : (np+2)/2;
    for (int i = 0; i < b.n; ++i) b.subnp[i] = np;
    for (int i = 0; i < b.n; ++i) b.offst[i] = 0;
    for (int min_np_pos = 0; min_np_pos < b.n; ++min_np_pos)
      recur(np, min_np, min_np_pos, good_bases, b, b.n-1, sa, basis, *mt);
    if ( ! good_bases.empty()) {
      max_good_np = std::max(max_good_np, min_np);
      ++good;
    }
    // We don't want to reduce order to improve other heuristics; instead, we'll
    // use the more general nodal subset basis if needed. So break at the
    // highest min_np having a good basis.
    if (good) break;
  }
  return max_good_np;
}

static void runall (const Int np,
                    const SearchAtom::Input::Basis basis = SearchAtom::Input::gll,
                    const MetricsTracker::Ptr& mt = nullptr) {
  if (np <= -1)
    for (int np : nps)
      run(np, mt, basis);
  else
    run(np, mt, basis);
}
} // namespace find_offset_nodal_subset_bases

// List of region's valid node supports. Use a vector of ints, where in an int,
// the 0/1 pattern gives the nodes.
struct ValidNodesList {
  typedef std::shared_ptr<ValidNodesList> Ptr;

  struct Iterator {
    Iterator (const std::vector<int>& nodes, const int& ntot, const int& nsub,
              const int& pin0)
      : nodes_(nodes), ntot_(ntot), nsub_(nsub), pin0_(pin0), idx_(-1)
    {}

    Iterator (const std::vector<int>& nodes, bool)
      : nodes_(nodes), ntot_(-1), nsub_(-1), pin0_(-1), idx_(nodes.size())
    {}

    Iterator operator++ () {
      int idx;
      for (idx = idx_+1; idx < static_cast<int>(nodes_.size()); ++idx)
        if (valid(nodes_[idx], pin0_))
          break;
      idx_ = idx;
      return *this;
    }

    bool operator== (const Iterator& it) const { return idx_ == it.idx_; }
    bool operator!= (const Iterator& it) const { return ! (*this == it); }

    // Convert a bit string to a list of element node indices.
    void get_nodes (int* nodes) {
      assert(idx_ < static_cast<int>(nodes_.size()));
      const auto nodemask = nodes_[idx_];
      int k = 0;
      for (int i = 0; i < ntot_; ++i)
        if (nodemask & (1 << i))
          nodes[k++] = i;
      assert(k == nsub_);
    }

  private:
    const std::vector<int>& nodes_;
    const int ntot_, nsub_, pin0_;
    int idx_;

    // Do the selected nodes include the two bounding the region?
    static bool valid (const int nodemask, const int pin0) {
      const int pinmask = 3 << pin0;
      return (nodemask & pinmask) == pinmask;
    }
  };

  ValidNodesList (const int ntot, const int nsub,
                  const bool symmetric = false) {
    assert(ntot <= islet::np_max);
    init(ntot, nsub, symmetric);
  }

  void init (const int ntot, const int nsub, const bool symmetric) {
    ntot_ = ntot;
    nsub_ = nsub;
    symmetric_ = symmetric;
    for (int i = 0; i < (1 << ntot); ++i)
      if (num1s(i) == nsub && ( ! symmetric_ || issymmetric(i, ntot)))
        nodes_.push_back(i);
    //pr(puf(ntot) pu(nsub) pu(nodes_.size()));
  }

  Iterator begin (const int& pin0) const {
    Iterator it(nodes_, ntot_, nsub_, pin0);
    ++it;
    return it;
  }

  Iterator end () const { return Iterator(nodes_, true); }

  static Int test () {
    if ( ! (issymmetric(1, 1) && issymmetric(3, 2) && !issymmetric(3, 3) &&
            issymmetric(0x65a6, 16) && !issymmetric(0x65a6, 20))) {
      std::cerr << "ValidNodesList::test FAILed.\n";
      return 1;
    }
    return 0;
  }

private:
  int ntot_, nsub_;
  bool symmetric_;
  std::vector<int> nodes_;

  // Number of 1 bits in n.
  static int num1s (int n) {
    int cnt = 0;
    while (n) {
      if (1 & n) ++cnt;
      n = n >> 1;
    }
    return cnt;
  }

  static bool issymmetric (const int n, const int nslot) {
    int nd = n, rev = 0;
    for (int i = 0; i < nslot; ++i) {
      rev = rev << 1;
      if (1 & nd) rev = (rev | 1);
      nd = nd >> 1;
    }
    return rev == n;
  }
};

namespace find_nodal_subset_bases {
struct NsbSearchAtom : public UserInterpMethod {
  struct Input {
    static const int np_max = islet::np_max;
    typedef char SInt;

    SInt np;
    SInt nodes[np_max-1][np_max];
    SInt subnp[np_max-1];
    Real maxeigampm1;
    Int ne, neigdx;
    bool quiet;

    Input () {
      np = -1;
      maxeigampm1 = 1e-13;
      ne = 1111;
      neigdx = ne;
      quiet = true;
    }
  };

  NsbSearchAtom (const bool mea_threaded = false)
    : max_eig_amp_(mea_threaded)
  {}

  Real run (const Input& in,
            bool& all_pve_wts, Real& wtr,
            MetricsTracker& mt, Real* metrics, Real& pum_metric) {
    in_ = in;
    Nodes nodes(in.np);
    for (Int i = 0; i < nodes.get_nh(); ++i)
      nodes.set(i, in.nodes[i], in.subnp[i]);
    assert(nodes.ok_to_eval());
    calc_xnodes_metrics(nodes, get_xnodes(), metrics);
    if ( ! mt.acceptable_metrics(nodes, get_xnodes(), metrics)) return 2;
    const Real pum_to_accept = mt.pum_to_accept(nodes, get_xnodes(), metrics);
    {
      Real wt[islet::np_max];
      calc_weights(nodes, get_xnodes(), wt);
      calc_wts_metrics(in.np, wt, all_pve_wts, wtr);
      if ( ! all_pve_wts) return 1;
    }
    Int ne, neigdx;
    ne = neigdx = 11;
    auto maxeigampm1 = max_eig_amp_.run(in.np, ne, neigdx, in.maxeigampm1,
                                        in.quiet, this);
    if (maxeigampm1 > in.maxeigampm1) return maxeigampm1;
    pum_metric = calc_pum_metric(*this, max_eig_amp_.is_threaded(),
                                 pum_to_accept);
    if ( ! mt.would_update(metrics, pum_metric)) return 2;
    return max_eig_amp_.run(in.np, in.ne, in.neigdx, in.maxeigampm1,
                            in.quiet, this);
  }

  void eval (const Real& x, Real* const v) override {
    eval(in_.np, in_.nodes, in_.subnp, x, v);
  }

  Int get_np () const override { return in_.np; }

  const Real* get_xnodes () const override {
    return islet::get_x_gll(in_.np);
  }

  static void eval (
    const Int& np, const Input::SInt nodes[][Input::np_max],
    const Input::SInt subnp[], const Real& x, Real* const v)
  {
    if (x > 0) {
      eval(np, nodes, subnp, -x, v);
      for (int i = 0; i < np/2; ++i)
        std::swap(v[i], v[np-i-1]);
      return;
    }
    const auto x_gll = islet::get_x_gll(np);
    Real xsub[Input::np_max], vsub[Input::np_max];
    for (Int i = 0; i < np-1; ++i) {
      if (i < np-2 && x > x_gll[i+1]) continue;
      if (subnp[i] == np) {
        eval_lagrange_poly(x_gll, np, x, v);
      } else {
#ifndef NDEBUG
        { // Subregion's nodes must be included for the basis to be
          // interpolatory.
          int fnd = 0;
          for (Int j = 0; j < subnp[i]; ++j) {
            const auto nij = nodes[i][j];
            if (nij == i || nij == i+1) ++fnd;
          }
          if (fnd != 2) {
            pr(puf(np) pu(i) pu(x));
            islet::prarr("nodes[i]", nodes[i], subnp[i]);
          }
          assert(fnd == 2);
        }
#endif
        for (Int j = 0; j < subnp[i]; ++j)
          xsub[j] = x_gll[(int) nodes[i][j]];
        // Lagrange polynomial basis.
        std::fill(v, v + np, 0);
        eval_lagrange_poly(xsub, subnp[i], x, vsub);
        for (Int j = 0; j < subnp[i]; ++j)
          v[(int) nodes[i][j]] = vsub[j];
      }
      break;
    }
  }

private:
  Input in_;
  MaxEigComputer max_eig_amp_;
};

static const int nps[] = {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16};

struct Basis {
  static const int np_max = NsbSearchAtom::Input::np_max;
  int np;
  std::array<int,np_max-1> subnp;
  std::array<std::array<int,np_max>,np_max-1> nodes;
};

struct Count {
  typedef std::shared_ptr<Count> Ptr;
  bool just_count;
  size_t value, total;
  Count () : just_count(true), value(0), total(0) {}
};

void randperm (Int* const p, const Int n) {
  using islet::urand;
  for (Int i = 0; i < n; ++i) p[i] = i;
  for (Int i = 0; i < 5*n; ++i) {
    const int j = urand()*n, k = urand()*n;
    std::swap(p[j], p[k]);
  }
}

void make_subnp_list_recur (const Int np, const Int np_min,
                            std::vector<char>& subnp_list,
                            char* subnp, const Int pos) {
  for (Int np_pos = np_min; np_pos <= np; ++np_pos) {
    subnp[pos] = np_pos;
    if (pos > 0)
      make_subnp_list_recur(np, np_min, subnp_list, subnp, pos-1);
    else
      for (Int j = 0; j < np/2; ++j)
        subnp_list.push_back(subnp[j]);
  }
}

void make_subnp_list (const Int np, const Int np_min,
                      std::vector<char>& subnp_list) {
  char subnp[islet::np_max];
  make_subnp_list_recur(np, np_min, subnp_list, subnp, np/2-1);
}

struct Restarter {
  typedef std::shared_ptr<Restarter> Ptr;
  MetricsTracker::Ptr mt;
  int np, min_ooa;
  size_t dont_eval_if_below, eval_count;

  Restarter (const MetricsTracker::Ptr& mt_, const int np_, const int min_ooa_,
             const int dont_eval_if_below_)
    : mt(mt_), np(np_), min_ooa(min_ooa_), dont_eval_if_below(dont_eval_if_below_),
      eval_count(0)
  {}

  bool write(std::string filename = "");
};

bool Restarter::write (std::string filename) {
  using islet::write;
  if (filename == "") {
    std::stringstream ss;
    ss << "NsbSearchAtomRestart_np" << np << ".dat";
    filename = ss.str();
  }
  std::ofstream os(filename.c_str(), std::ofstream::binary);
  assert(mt);
  return (write(os, np) &&
          write(os, min_ooa) &&
          write(os, eval_count) && // new dont_eval_if_below value
          write(os, eval_count) &&
          mt->write(os));
}

static Restarter::Ptr read_restart (const int np) {
  using islet::read;
  std::stringstream ss;
  ss << "NsbSearchAtomRestart_np" << np << ".dat";
  std::ifstream is(ss.str().c_str(), std::ofstream::binary);
  if ( ! is.is_open()) return nullptr;
  const auto mt = std::make_shared<MetricsTracker>(np);
  const auto r = std::make_shared<Restarter>(mt, np, 0, 0);
  const bool ok = (read(is, r->np) && r->np == np &&
                   read(is, r->min_ooa) &&
                   read(is, r->dont_eval_if_below) &&
                   read(is, r->eval_count) &&
                   r->mt->read(is));
  if ( ! ok) return nullptr;
  return r;
}

void eval (std::vector<NsbSearchAtom>& esa,
           const std::vector<NsbSearchAtom::Input>& input_list,
           const Int ninput, MetricsTracker& mt, const bool show_progress,
           const Count::Ptr& count) {
  if (count) count->value += ninput;
  if (count && count->just_count) return;
  if (count)
    printf("NsbSearchAtom::eval %ld/%ld (%5.1f%%)\n",
           count->value, count->total, 100*Real(count->value)/count->total);
  std::vector<Int> perm(ninput);
  randperm(perm.data(), ninput);
  Real progress = 0;
  const auto run1 = [&] (const Int ili) {
    {
      const Real pdelta = 0.05;
      const Real p = Real(ili)/ninput;
      if (show_progress && p >= progress + pdelta) {
#       pragma omp critical (NsbSearchAtom_progress)
        {
          const Real p = Real(ili)/ninput;
          if (p >= progress + pdelta) {
            printf("progress: %5.1f%% (%8d)\n", 100*p, ili);
            progress = p;
          }
        }
      }
    }
    const auto& in = input_list[perm[ili]];
    auto& my_esa = esa[omp_get_thread_num()];
    bool all_pve_wts = false;
    Real wtr, metrics[3], pum_metric;
    const auto maxeigampm1 = my_esa.run(in, all_pve_wts, wtr, mt, metrics,
                                        pum_metric);
    if (maxeigampm1 > in.maxeigampm1) return;
    {
      Basis b;
      b.np = in.np;
      for (int i = 0; i < b.np-1; ++i)
        b.subnp[i] = in.subnp[i];
      for (int i = 0; i < b.np-1; ++i)
        for (int j = 0; j < b.subnp[i]; ++j)
          b.nodes[i][j] = in.nodes[i][j];
#     pragma omp critical (NsbSearchAtom_eval)
      {
        mt.update(metrics, pum_metric);
        const int n = in.np/2;
        printf("meam1 %9.2e w>0 %d wtr %8.2e npm %8.2e %8.2e %8.2e pum %9.2e | ",
               maxeigampm1, all_pve_wts, wtr, metrics[0], metrics[1], metrics[2],
               pum_metric);
        printf("np %2d subnp", in.np);
        for (int i = 0; i < n; ++i) printf(" %d", in.subnp[i]);
        printf(" nodes");
        for (int i = 0; i < n; ++i) {
          printf(" |");
          for (int j = 0; j < in.subnp[i]; ++j)
            printf(" %d", in.nodes[i][j]);
        }
        printf("\n");
      }
    }
  };

  if (esa.size() > 1) {
#   pragma omp parallel for schedule(dynamic,1)
    for (Int ili = 0; ili < ninput; ++ili)
      run1(ili);
  } else {
    for (Int ili = 0; ili < ninput; ++ili)
      run1(ili);
  }
}

void recur (const int np, const std::vector<ValidNodesList::Ptr>& vnls,
            const std::vector<ValidNodesList::Ptr>& vnls_mid_reg,
            Basis b, const int pos, std::vector<NsbSearchAtom>& esa,
            std::vector<NsbSearchAtom::Input>& input_list,
            Int& input_list_pos, MetricsTracker& mt, const Count::Ptr& count,
            Restarter& restarter) {
  if (pos == -1) {
    NsbSearchAtom::Input& in = input_list[input_list_pos];
    in.np = np;
    in.maxeigampm1 = 1e-13;
    for (int i = 0; i < b.np-1; ++i)
      in.subnp[i] = b.subnp[i];
    for (int i = 0; i < b.np-1; ++i)
      for (int j = 0; j < b.subnp[i]; ++j)
        in.nodes[i][j] = b.nodes[i][j];
    ++input_list_pos;
    if ((size_t) input_list_pos == input_list.size()) {
      if (restarter.eval_count >= restarter.dont_eval_if_below) {
        // Run a bunch of analyses in parallel.
        eval(esa, input_list, input_list_pos, mt, false, count);
        if ( ! (count && count->just_count)) {
          if ( ! count)
            printf("restart eval_count %ld\n", restarter.eval_count);
          ++restarter.eval_count;
          restarter.write();
        }
      } else {
        if (count) count->value += input_list_pos;
        if ( ! (count && count->just_count)) ++restarter.eval_count;
      }
      input_list_pos = 0;
    }
    return;
  }
  // Set up a basis.
  const auto& vnl = (np % 2 == 0 && pos == np/2-1) ?
    vnls_mid_reg[b.subnp[pos]] : vnls[b.subnp[pos]];
  for (auto it = vnl->begin(pos); it != vnl->end(); ++it) {
    it.get_nodes(b.nodes[pos].data());
    recur(np, vnls, vnls_mid_reg, b, pos-1, esa, input_list, input_list_pos,
          mt, count, restarter);
  }
}

static Int run (const int np, int min_ooa = -1,
                MetricsTracker::Ptr mt = nullptr,
                const Count::Ptr count = nullptr,
                const int dont_eval_if_below = 0) {
  assert(np <= NsbSearchAtom::Input::np_max);
  if ( ! mt) mt = std::make_shared<MetricsTracker>(np);
  const bool thread_toplevel = np >= 9;
  std::vector<NsbSearchAtom> esa;
  if (thread_toplevel)
    esa.resize(omp_get_max_threads());
  else
    esa.emplace_back(true);
  std::vector<NsbSearchAtom::Input> input_list(thread_toplevel ?
                                               1 << 22 :
                                               1 << 12);
  Int input_list_pos = 0;
  const bool just_count = count && count->just_count;
  if ( ! just_count) printf("np %2d\n", np);
  if (min_ooa <= 0) {
    if (np == 5)      min_ooa = 2;
    else if (np <= 6) min_ooa = np-2;
    else if (np <= 9) min_ooa = np-3;
    else              min_ooa = np-4;
  }
  std::vector<ValidNodesList::Ptr> vnls(np+1), vnls_mid_reg(np+1);
  for (int min_np = np; min_np >= min_ooa+1; --min_np)
    vnls[min_np] = std::make_shared<ValidNodesList>(np, min_np);
  for (int min_np = np; min_np >= min_ooa+1; --min_np)
    vnls_mid_reg[min_np] = std::make_shared<ValidNodesList>(np, min_np, true);
  std::vector<char> subnp_list;
  make_subnp_list(np, min_ooa+1, subnp_list);
  Restarter restarter(mt, np, min_ooa, dont_eval_if_below);
  const Int sz = np/2;
  const Int n = subnp_list.size()/sz;
  for (Int isubnp = 0; isubnp < n; ++isubnp) {
    const char* subnp = &subnp_list[isubnp*sz];
    Basis b;
    b.np = np;
    for (Int j = 0; j < sz; ++j) {
      b.subnp[j] = subnp[j];
      b.subnp[np-2-j] = subnp[j];
    }
    recur(np, vnls, vnls_mid_reg, b, b.np/2 - 1, esa,
          input_list, input_list_pos, *mt, count, restarter);
  }
  if (input_list_pos > 0)
    eval(esa, input_list, input_list_pos, *mt, true, count);
  return 0;
}

static void run () {
  for (int np : nps) run(np);
}
} // namespace find_nodal_subset_bases

static void find_nodal_subset_bases_given_mt (
  const int np, const MetricsTracker::Ptr& mt, const int min_ooa = -1,
  const int eval_count = 0)
{
  if (np > 10) {
    find_nodal_subset_bases::run(np, min_ooa, mt, nullptr, eval_count);
    return;
  }
  const auto count = std::make_shared<find_nodal_subset_bases::Count>();
  find_nodal_subset_bases::run(np, min_ooa, mt, count);
  printf("count %ld\n", count->value);
  count->total = count->value;
  count->value = 0;
  count->just_count = false;
  find_nodal_subset_bases::run(np, min_ooa, mt, count, eval_count);
}

static void find_nodal_given_best_offset_nodal (
  const int np, const bool restart_if_available = true)
{
  find_nodal_subset_bases::Restarter::Ptr restarter;
  if (restart_if_available)
    restarter = find_nodal_subset_bases::read_restart(np);
  if (restarter) {
    find_nodal_subset_bases_given_mt(np, restarter->mt, restarter->min_ooa,
                                     restarter->eval_count);
  } else {
    const auto mt = std::make_shared<MetricsTracker>(np);
    const Int max_good_np = find_offset_nodal_subset_bases::run(np, mt);
    const Int min_ooa = max_good_np-1;
    mt->set_pum_max(mt->get_pum_min());
    find_nodal_subset_bases_given_mt(np, mt, min_ooa);
  }
}

static void run_general_unittests () {
  int nerr = 0;
  {
    for (int np = 2; np <= 7; ++np) {
      const auto x = islet::get_x_gll(np);
      const auto w = islet::get_w_gll(np);
      Real sum = 0;
      for (int j = 0; j < np; ++j) sum += w[j];
      if (islet::reldif(2, sum) >= 1e-14) ++nerr;
      for (int j = 0; j < np/2; ++j)
        if (w[j] != w[np-j-1]) ++nerr;
      for (int j = 0; j < np/2; ++j)
        if (x[j] != -x[np-j-1]) ++nerr;
      for (int j = 0; j < np-1; ++j)
        if (x[j+1] < x[j]) ++nerr;
    }
  }
  nerr += MaxEigComputer::unittest();
  nerr += SearchAtom::unittest();
  nerr += ValidNodesList::test();
  std::cout << (nerr ? "FAIL" : "PASS") << " unit test\n";
}

namespace find_natural_bases {
// Sweep over node locations and check the natural basis for positive weights
// and, if positive, stability.

struct NaturalInterp : public UserInterpMethod {
  Int np = 0;
  Real xnodes[islet::np_max] = {0};
  
  virtual void eval (const Real& x, Real* const v) override {
    eval_lagrange_poly(xnodes, np, x, v);
  }

  virtual const Real* get_xnodes() const override { return xnodes; }

  virtual Int get_np() const override { return np; }
};

struct NaturalXnodes {
  NaturalXnodes (const Int inp, const Int insample)
    : np(inp), tol(1e-13), ne_max(500), ndx_max(500), nsample(insample),
      nslot(np/2-1), mt(np), verbose(false)
  {
    nodes.init(np, true);
    std::vector<Int> subset(np);
    for (int i = 0; i < np; ++i) subset[i] = i;
    for (int i = 0; i < nodes.get_nh(); ++i) nodes.set(i, subset);
    assert(nodes.ok_to_eval());

    im.np = np;
    im.xnodes[0] = -1;
    im.xnodes[np-1] =  1;
    if (np % 2 == 1) im.xnodes[np/2] = 0;
  }

  void set_verbose (const bool b) { verbose = b; }

  void run () {
    recur(0, 0);
    print_stats();
  }

  void print_stats () {
    printf("neval %d npassfilt %d npumeval %d n_mt_noaccept %d\n",
           neval, npassfilt, npumeval, n_mt_noaccept);
  }

private:
  Real tol;
  Int np, ne_max, ndx_max, nsample, nslot;
  Nodes nodes;
  MetricsTracker mt;
  NaturalInterp im;
  MaxEigComputer max_eig_amp;
  bool verbose;
  Int neval = 0, npassfilt = 0, npumeval = 0, n_mt_noaccept = 0;

  void recur (const Int slot, const Real x_base) {
    if (slot == nslot) {
      ++neval;

      assert_stuff();

      const auto* xnodes = im.get_xnodes();
      Real xnodes_metric[3];
      calc_xnodes_metrics(nodes, xnodes, xnodes_metric);
      if ( ! mt.acceptable_metrics(nodes, xnodes, xnodes_metric)) {
        ++n_mt_noaccept;
        return;
      }
      bool all_pve_wts = false;
      Real wtr = 0; {
        Real wt[islet::np_max];
        calc_weights(nodes, xnodes, wt);
        calc_wts_metrics(np, wt, all_pve_wts, wtr);
      }
      if ( ! all_pve_wts) return;
      ++npassfilt;

      const Real maxeigampm1 = max_eig_amp.run(np, ne_max, ndx_max, tol, true,
                                               &im);
      const auto good = maxeigampm1 <= tol;
      Real pum_metric = 1;
      if (good) {
        ++npumeval;
        const Real pum_to_accept = mt.pum_to_accept(nodes, xnodes, xnodes_metric);
        pum_metric = calc_pum_metric(im, pum_to_accept);
        mt.update(xnodes_metric, pum_metric);
      }

      if (good || verbose) {
        printf(" (%1.3e, %1.3e, %1.3e) %s",
               maxeigampm1, pum_metric, wtr, good ? "good" : "");
        printf(" x");
        for (int i = 0; i < np; ++i)
          printf(" %22.15e", xnodes[i]);
        printf("\n");
      }
    } else {
      const int n = std::ceil(nsample*(1 - x_base)) + 1;
      assert(n >= 2);
      for (int i = 1; i < n; ++i) {
        const Real a = Real(i)/n;
        const Real x = (1-a)*x_base + a;
        im.xnodes[(np+1)/2+slot] = x;
        im.xnodes[np/2-1-slot] = -x;
        recur(slot+1, x);
      }
    }    
  }

  void assert_stuff () {
#ifndef NDEBUG
    const auto* x = im.get_xnodes();
    for (int i = 1; i < np; ++i)
      assert(x[i] > x[i-1]);
    for (int i = 0; i < np; ++i)
      assert(x[i] == -x[np-1-i]);
#endif
  }
};

void run (const Int np, const Int nsample) {
  NaturalXnodes nx(np, nsample);
  nx.run();
}
} // namespace find_natural_bases

struct Command {
  enum Enum { unittest, findoffsetnodal, findnodal,
              finduniform, findlegendre, findcheb,
              findnodal_given_bestosn, findnatural};
  static Enum convert (const std::string& s) {
    if (s == "unittest") return unittest;
    if (s == "findoffsetnodal") return findoffsetnodal;
    if (s == "findnodal") return findnodal;
    if (s == "findnodal_given_bestosn") return findnodal_given_bestosn;
    if (s == "finduniform") return finduniform;
    if (s == "findlegendre") return findlegendre;
    if (s == "findcheb") return findcheb;
    if (s == "findnatural") return findnatural;
    throw std::logic_error("Not a command.");
  }
};

int main (int argc, char** argv) {
  if (argc < 2) {
    std::cerr << argv[0] << " <command> <options>\n";
    return -1;
  }

  using NosbBasis = SearchAtom::Input;

  const auto command = Command::convert(argv[1]);
  bool unittest = false;
  Int np = -1;
  if (argc > 2) np = std::atoi(argv[2]);
  switch (command) {
  case Command::unittest: {
    unittest = true;
  } break;
  case Command::findoffsetnodal: {
    // Pretty efficient search for stable offset nodal bases.
    find_offset_nodal_subset_bases::runall(np);
  } break;
  case Command::findnodal: {
    // Search for stable nodal subset bases.
    const auto restarter = find_nodal_subset_bases::read_restart(np);
    if (restarter)
      find_nodal_subset_bases_given_mt(np, restarter->mt, restarter->min_ooa,
                                       restarter->eval_count);
    else
      find_nodal_subset_bases_given_mt(np, nullptr);
  } break;
  case Command::findnodal_given_bestosn: {
    if (np == -1) {
      std::cerr << argv[0] << " findnodal_given_bestosn np\n";
      return -1;      
    }
    find_nodal_given_best_offset_nodal(np);
  } break;
  case Command::finduniform:
    find_offset_nodal_subset_bases::runall(np, NosbBasis::uniform); break;
  case Command::findlegendre:
    find_offset_nodal_subset_bases::runall(np, NosbBasis::legendre); break;
  case Command::findcheb:
    find_offset_nodal_subset_bases::runall(np, NosbBasis::cheb); break;
  case Command::findnatural: {
    int nsample = 111;
    if (argc > 3) nsample = std::atoi(argv[3]);
    if (argc >= 1)
      find_natural_bases::run(np, nsample);
  } break;
  default:
    throw std::logic_error("Not a command.");
  }
  if (unittest) run_general_unittests();
}
