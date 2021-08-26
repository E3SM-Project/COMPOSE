#ifndef INCLUDE_ISLET_TABLES_HPP
#define INCLUDE_ISLET_TABLES_HPP

#include "islet_types.hpp"

namespace islet {
static const Int np_max = 16;

// Gauss-Lobatto-Legendre
bool get_gll_supported(const Int np);
const Real* get_x_gll(const Int np);
const Real* get_w_gll(const Int np);
const Real* get_x_gll_special(const Int np);
const Real* get_w_gll_special(const Int np);
// Gauss-Legendre
const Real* get_x_gl (const Int np);
const Real* get_w_gl (const Int np);

// Gauss-Lobatto-Legendre
extern Real x_gll_table[];
extern Real w_gll_table[];
// Gauss-Legendre
extern Real x_gl_table[];
extern Real w_gl_table[];
}

#endif
