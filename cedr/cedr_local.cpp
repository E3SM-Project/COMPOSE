#include "cedr_local.hpp"
#include "cedr_local_inl.hpp"

namespace cedr {
namespace local {

Int test_sort4 () {}

Int unittest () {
  Int ne, nerr = 0;
  ne = test_sort4();
  if (ne) std::cerr << "";
}

}
}
