#ifndef INCLUDE_MPAS_LAGRANGIAN_LEVELS_HPP
#define INCLUDE_MPAS_LAGRANGIAN_LEVELS_HPP

#include "mpas_mesh.hpp"
#include "mpas_trajectory.hpp"
#include "mpas_ci_util.hpp"
#include "mpas_comm.hpp"

#include <memory>

struct LagrangianLevels {
  typedef std::shared_ptr<LagrangianLevels> Ptr;

  struct Params {
    bool cell_integrated, bfb;
    Params () : cell_integrated(true), bfb(true) {}
  };

  LagrangianLevels(const MpasMesh::ConstPtr& mm, const siqp::sh::Mesh::ConstPtr& m,
                   const Comm::Ptr& he, const Trajectory::ConstPtr& traj,
                   // If traj->has_in_cell_data(), otm is unused.
                   const OctreeMgr::ConstPtr& otm, const Params p = Params());

  // wrk must have size
  //     order_invariant_atomic_accum_size_nreal*nilvl*ncell.
  void map_zeta_to_eulerian_cell_centers(
    const Acr3& p_advect_vtx, const Ar3& zeta_vtx, const Ar1& wrk,
    const Ar2& zeta_lag_cc);

  void remap_vertically(const Acr2& zeta_eul, const Ar2& zeta_lag, const Ar2& rho,
                        const Ar1& rwrk, const int nq, const Ar3& Q);

  static bool unittest(
    const MpasMesh::Ptr& mm, const siqp::sh::Mesh::ConstPtr& m, const Comm::Ptr& he,
    const bool verbose = false);

private:
  struct Impl;
  std::shared_ptr<Impl> impl;
};

#endif
