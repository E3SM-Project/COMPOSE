#include "mpas_trajectory_impl.hpp"

#include <map>

static Real get_dist2 (const Real* const a, const Real* const b) {
  Vec amb; geo::axpbyz(1, a, -1, b, amb);
  return geo::norm2(amb);
}

// A cove is an inlet in this rank's cells (owned + shared) in which a cell is
// closer to a point p than all its neighbors, without actually containing p. We
// need to escape the cove to continue search for the cell that contains p.
int Trajectory::Impl
::escape_cove (const int ci, Real min_dist2, const Real* const p) {
  const auto& dist = he->cdhalo();
  const auto& ccs = mm->g.my_cellCentroid;
  const int nnbr = c2ccnt[ci];
  const int* ci_c2c = c2c(ci);
  int min_cj = ci;
  for (int j = 0; j < nnbr; ++j) {
    const int cj = ci_c2c[j];
    assert(cj >= 0);
    if (dist(cj) == 1) continue;
    Vec dp;
    geo::axpbyz(1, p, -1, ccs(cj), dp);
    const Real d2 = geo::norm2(dp);
    if (d2 < min_dist2) {
      min_dist2 = d2;
      min_cj = cj;
    }
    if (min_cj != ci) break;
  }
  assert(min_cj >= 0);
  return min_cj;
}

int Trajectory::Impl
::find_cell_given_starting_cell (const int ci0, const Real* const p, int* pcnt) {
  const auto& t = mm->t;
  const auto& ccs = mm->g.my_cellCentroid;
  const auto& dist = he->cdhalo();
  int ci = ci0, cnt = 0;
  for (;;) {
    ++cnt;
    Real min_dist2 = get_dist2(ccs(ci), p);
    int min_ci = ci;
    const auto ne = t.c2ne[ci];
    const int* nbrs = c2c(ci);
    for (int j = 0; j < ne; ++j) {
      const int cj = nbrs[j];
      if (dist(cj) == 1) continue;
      const Real d2 = get_dist2(ccs(cj), p);
      if (d2 < min_dist2) {
        min_dist2 = d2;
        min_ci = cj;
      }
    }
    if (min_ci == ci) {
      if (min_dist2 < cc_min_nbr_dist2(ci)) break;
      // Otherwise, we're currently stuck in a cove. Broaden the search to get
      // unstuck.
      min_ci = escape_cove(ci, min_dist2, p);
      if (min_ci == ci) {
        // Give up; we can't find a closer cell.
        break;
      }
    }
    ci = min_ci;
  }
  if (pcnt) *pcnt = cnt;
  return ci;
}

static void
make_point (Trajectory::Impl& ti, const int ci, const int place, Vec p) {
  const auto& t = ti.mm->t;
  const auto& g = ti.mm->g;
  const auto& v = ti.mm->v;
  switch (place) {
  case 0: { // cell ctr
    geo::copy(p, g.my_cellCentroid(ci));
  } break;
  case 1: { // random pt in cell
    const auto ne = t.c2ne[ci];
    for (int d = 0; d < 3; ++d) p[d] = 0;
    Real a[max_vtx], den = 0;
    for (int i = 0; i < ne; ++i) a[i] = urand();
    for (int i = 0; i < ne; ++i) den += a[i];
    for (int i = 0; i < ne; ++i) a[i] /= den;
    for (int i = 0; i < ne; ++i) {
      const int vi = t.get_c2v(ci,i);
      p[0] += a[i]*v.x[vi];
      p[1] += a[i]*v.y[vi];
      p[2] += a[i]*v.z[vi];
    }
    geo::normalize(p);
  } break;
  }  
}

static int
run_search_test (Trajectory::Impl& ti, const int cib, const int cie,
                 const int place, Real& dist, int& cnt) {
  Vec p;
  make_point(ti, cie, place, p);
  const auto& ccs = ti.mm->g.my_cellCentroid;
  const Real d2b = std::sqrt(get_dist2(ccs(cib), p));
  dist = d2b;
  const auto ci_search = ti.find_cell_given_starting_cell(cib, p, &cnt);
  assert(ci_search >= 0);
  assert(ci_search < ti.mm->t.ncell);
  if (ci_search != cie) {
    const Real d2e = std::sqrt(get_dist2(ccs(cie), p));
    const Real d2s = std::sqrt(get_dist2(ccs(ci_search), p));
    pr("run_search_test FAIL" pu(cib) pu(cie) pu(ci_search) pu(place) pu(d2b)
       pu(d2e) pu(d2s));
    return 1;
  }
  return 0;
}

struct HistBin {
  int n;
  Real lo, hi;
  HistBin (const Real dist) : n(1), lo(dist), hi(dist) {}
};

static Int get_glbl_num_cells (const MpasMesh& mm, Comm& he) {
  Int gncell;
  he.all_reduce(Comm::SUM, 1, &mm.ncown, &gncell);
  return gncell;
}

int Trajectory::test_search (const MpasMesh::ConstPtr& mm, const Comm::Ptr& he,
                             const bool print_stats) {
  int nerr = 0;
  const auto vf = vf_from_short('n');
  const auto traj = std::make_shared<Trajectory>(
    mm, he, vf, Trajectory::Method::taylor2, 6);
  const auto& he_dist = he->cdhalo();
  const auto ncown = mm->ncown;
  std::vector<Real> dists;
  std::vector<int> cnts;
  for (int trial = 0; trial < 100; ++trial) {
    for (const int place : {0, 1}) {
      int cib, cie;
      do {
        cib = int(urand()*ncown);
        cie = int(urand()*ncown);
      } while (he_dist(cib) == 1 || he_dist(cie) == 1);
      int cnt;
      Real dist;
      nerr += run_search_test(*traj->impl, cib, cie, place, dist, cnt);
      if (print_stats) {
        cnts.push_back(cnt);
        dists.push_back(dist);
      }
    }
  }
  // 4 pi r^2 = ncell => 2 pi r =
  const auto gncell = get_glbl_num_cells(*mm, *he);
  const auto circum = 2*M_PI*std::sqrt(gncell/(4*M_PI));
  if (print_stats) {
    std::map<int, HistBin> cnts_hist;
    size_t ncnt = cnts.size();
    for (size_t i = 0; i < ncnt; ++i) {
      const auto cnt = cnts[i];
      auto it = cnts_hist.find(cnt);
      auto d = dists[i];
      // Convert from Euclidean distance to spherical distance / pi.
      d = (2*std::asin(d/2))/M_PI;
      if (it == cnts_hist.end())
        cnts_hist.emplace(cnt, HistBin(d));
      else {
        auto& hb = it->second;
        hb.n++;
        hb.lo = std::min(hb.lo, d);
        hb.hi = std::max(hb.hi, d);
      }
    }
    printf("circumference in cells: %5.1f\n", circum);
    for (const auto& kv : cnts_hist) {
      const auto cnt = kv.first;
      const auto& hb = kv.second;
      printf("%3d: %5.1f%% %5.3f %5.3f\n",
             cnt, 100*Real(hb.n)/ncnt, hb.lo, hb.hi);
    }
  }
  return nerr;
}
