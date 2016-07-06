#include "slmm_gallery.hpp"

namespace slmm {
namespace gallery {

const char* InitialCondition::inputs[] =
  {"xyztrig", "gaussianhills", "cosinebells", "slottedcylinders",
   "correlatedcosinebells"};

const char* WindFieldType::inputs[] =
  {"dcmip1d3ll", "nondivergent", "divergent", "rotate", "nondivergenthack"};

} // namespace gallery
} // namespace slmm
