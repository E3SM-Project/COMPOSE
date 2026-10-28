#ifndef INCLUDE_MPAS_ZETA_HPP
#define INCLUDE_MPAS_ZETA_HPP

#include "mpas_mesh.hpp"
#include "mpas_gallery.hpp"
#include "mpas_ci_util.hpp"

struct ZetaRep {
  typedef std::shared_ptr<const ZetaRep> ConstPtr;

  Acr2 eul_cell;            // (li,ci)
  Ar2  eul_vtx_univalued;   // (li,vi)
  Ar3  eul_vtx_multivalued; // (vtx,li,ci)

  ZetaRep(const MpasMesh::ConstPtr& mm, const siqp::sh::Mesh& m,
          const Comm::Ptr& he);
};

#endif
