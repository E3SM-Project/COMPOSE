#include "mpas_debug.hpp"

namespace dbg {
void write_array (const Array1D<const int>* gids, const Array3D<const double>& d,
                  const int d0, const int d1, const int d2beg, const int d2end,
                  const std::string& filename) {
  write_array_impl(gids, d, d0, d1, d2beg, d2end, filename);
}
} // namespace dbg
