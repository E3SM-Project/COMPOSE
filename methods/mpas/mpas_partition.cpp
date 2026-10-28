#include "mpas_partition.hpp"
#include "mpas_util.hpp"

#include <fstream>
#include <sstream>

std::vector<Int> read_partition (const std::string& part_file) {
  std::ifstream fp;
  fp.open(part_file, std::ifstream::in);
  std::vector<Int> c2r;
  while (fp.good()) {
    int i;
    fp >> i;
    if (fp.rdstate() & std::fstream::failbit) break;
    c2r.push_back(i);
  }
  return c2r;
}

std::vector<Int> create_serial_partition (const Int ncell) {
  return std::vector<Int>(ncell, 0);
}

// r2cs[r][i] is the ith cell rank r owns, in LID order.
std::vector<std::vector<Int> > c2r_to_r2cs (const std::vector<Int>& c2r) {
  Int rmax = 0;
  for (const auto& r : c2r) rmax = std::max(rmax, r);
  std::vector<std::vector<Int> > r2cs(rmax+1);
  const Int nc = c2r.size();
  for (Int c = 0; c < nc; ++c) r2cs[c2r[c]].push_back(c);
  return r2cs;
}

CellGraph::ConstPtr read_cell_graph (const std::string& grid_file) {
  CellGraph::Ptr pcg = std::make_shared<CellGraph>();
  auto& cg = *pcg;
  std::ifstream fp;
  fp.open(grid_file, std::ifstream::in);
  if ( ! fp.good()) return pcg;
  Int ncell, nedge; {
    std::string line;
    std::getline(fp, line);
    std::stringstream ss(line);
    if ( ! (ss >> ncell >> nedge)) return pcg;
    nedge *= 2;
  }
  cg.cells.resize(ncell+1);
  cg.nbrs.resize(nedge);
  cg.cells[0] = 0;
  for (Int ci = 0; ci < ncell; ++ci) {
    std::string line;
    std::getline(fp, line);
    std::stringstream ss(line);
    cg.cells[ci+1] = cg.cells[ci];
    for (;;) {
      Int c;
      if ( ! (ss >> c)) break;
      cg.nbrs[cg.cells[ci+1]++] = c-1;
    }
  }
#ifndef NDEBUG
  assert(cg.cells[ncell] == nedge);
  for (Int ci = 0; ci < ncell; ++ci) {
    const Int n = cg.cells[ci+1] - cg.cells[ci];
    assert(n >= 5);
    assert(n <= 8);
  }
  for (Int ei = 0; ei < nedge; ++ei) {
    assert(cg.nbrs[ei] >= 0);
    assert(cg.nbrs[ei] < ncell);
  }
#endif
  return pcg;
}
