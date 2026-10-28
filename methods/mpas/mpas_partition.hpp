#ifndef INCLUDE_MPAS_PARTITION_HPP
#define INCLUDE_MPAS_PARTITION_HPP

#include "mpas_mesh.hpp"

// Read the grid partition file from gpmetis. Assign cell i to rank c2r[i].
std::vector<Int> read_partition(const std::string& part_file);

std::vector<Int> create_serial_partition(const Int ncell);

// r2cs[r][i] is the ith cell rank r owns, in LID order.
std::vector<std::vector<Int> > c2r_to_r2cs(const std::vector<Int>& c2r);

// Indexing is base-0. Cell i's neighbors are nbrs(cells(i):cells(i+1)-1).
struct CellGraph {
  typedef std::shared_ptr<CellGraph> Ptr;
  typedef std::shared_ptr<const CellGraph> ConstPtr;
  std::vector<Int> cells;
  std::vector<Int> nbrs;
  int ncell () const { return cells.size() - 1; }
};

CellGraph::ConstPtr read_cell_graph(const std::string& grid_file);

#endif
