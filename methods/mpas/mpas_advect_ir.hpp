#ifndef INCLUDE_MPAS_ADVECT_IR_HPP
#define INCLUDE_MPAS_ADVECT_IR_HPP

#include "mpas_mesh.hpp"
#include "mpas_comm.hpp"
#include "mpas_trajectory.hpp"

#include <string>

namespace ir {

struct Params {
  VelocityField::Ptr vf;
  int ntracer, order, tq_order, nstep, nsubstep, traj_nhalo, limit, cedr;
  bool perturb_dep_pts; // add noise to departure points to test robustness
  bool lag_lvls_cell_int; // cell-integration or interp in zeta remap
  // Where there are two methods, one BFB w.r.t. MPI/thread and one not, use the
  // BFB one.
  bool bfb;
  Trajectory::Method::Enum advect_method;

  Params ()
    : ntracer(-1), order(2), tq_order(-1), nstep(-1), nsubstep(-1),
      traj_nhalo(-1), limit(0), cedr(1), perturb_dep_pts(0),
      lag_lvls_cell_int(true), bfb(true), advect_method(Trajectory::Method::order1)
  {}
};

struct RemapFormOp {
  typedef std::shared_ptr<RemapFormOp> Ptr;
  
  RemapFormOp(const MpasMesh::ConstPtr& mm, const siqp::sh::Mesh::ConstPtr& m,
              const Comm::Ptr& he, const Params& p);

  // If is_3d, rhos and Qs are modified on output.
  void integrate(const Real ts, const Real dt, const Ar2& rhos, const Ar2& rhot,
                 const Ar3& Qs, const Ar3& Qt);

  bool unittest(Real* data = nullptr);
  void advect_mesh(const Real ts, const Real dt);
  void fill_patches(std::vector<Real>& xyzs_orig, std::vector<Real>& xyzs_advected);

  MpasMesh::ConstPtr get_mpas_mesh() const;
  siqp::sh::Mesh::ConstPtr get_siqp_mesh() const;
  Comm::Ptr get_comm() const;
  const Params& get_params() const;

  static RemapFormOp::Ptr
  create(const MpasMesh::ConstPtr& mm, const siqp::sh::Mesh::ConstPtr& m,
         const Comm::Ptr& he, const Params& p);

private:
  struct Impl;
  std::shared_ptr<Impl> impl;
};

} // namespace ir

#endif
