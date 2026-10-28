#ifndef INCLUDE_MPAS_CHECK_HPP
#define INCLUDE_MPAS_CHECK_HPP

#include "mpas_mesh.hpp"

void check(const MpasMesh& m);
void get_angleEdge_error(const Topo& t, const CellCtrs& c, const EdgeCtrs& e,
                         const Vtxs& v, const MpasGeometry& g, Real* err);

#endif
