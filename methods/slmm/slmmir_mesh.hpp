#ifndef INCLUDE_SLMMIR_MESH_HPP
#define INCLUDE_SLMMIR_MESH_HPP

#include <memory>

#include "slmm_array.hpp"
#include "slmm_basis.hpp"
using namespace slmm;

struct Mesh {
  typedef std::shared_ptr<Mesh> Ptr;

  Int np, nonuni, tq_order;
  Basis::Ptr basis;
  AVec3s geo_p, geo_nml, cgll_p;
  AIdxs geo_c2n, geo_c2nml, cgll_c2n, dgll_c2n, cgll_io_c2n;
  AIdxArray dglln2cglln;
};

#endif
