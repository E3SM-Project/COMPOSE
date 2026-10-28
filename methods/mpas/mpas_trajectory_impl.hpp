#ifndef INCLUDE_MPAS_TRAJECTORY_IMPL_HPP
#define INCLUDE_MPAS_TRAJECTORY_IMPL_HPP

#include "mpas_trajectory.hpp"
#include "mpas_geometry.hpp"
#include "mpas_comm_util.hpp"
#include "mpas_ci_util.hpp"
#include "mpas_zeta.hpp"

struct Trajectory::Impl {
  static constexpr int e2v_function = 1;
  // Set at init.
  Method::Enum method;
  VelocityField::Ptr vf;
  MpasMesh::ConstPtr mm;
  Comm::Ptr he;
  Real method_exact_rel_tol;
  bool is_3d;
  int nlvl, nsubstep, nhalo, vdegree, c2cslotsz, v2cslotsz;
  Ai1 c2ccnt, v2ccnt;
  Ai2 c2c, v2c, e2v_stencils;
  Ar1 angle_edge, cc_min_nbr_dist2;
  Ar2 ecs, e2v_coefs;
  CellFromVtxsShapeFunction cvsf;
  // Mod'ed during time stepping.
  Array2D<CellLvl> vtx_in_cell_lvl;
  Ac2 edge_mask, vtx_mask, cell_mask;
  Ar1 initial_step;

  size_t get_rwrk_size() const;

  void advect_vtxs(const Acr2& p, const Real ts, const Real dt, const Ar1& rwrk,
                   const Ar3& p_advect);
  bool unittest(const Acr2& vtxs, Real* data = nullptr);
  static bool convtest(const MpasMesh::ConstPtr& m, const Comm::Ptr& he,
                       Real* data = nullptr);

  void init(const MpasMesh::ConstPtr& m, const Comm::Ptr& he, VelocityField::Ptr vf,
            const Method::Enum method, const int nsubstep, const int nhalo);
  void setup_advect_vtxs();

  void masks_init_from_own(const char maskval, const bool cell_halo);
  void masks_init_from_cells(const Array2D<const CellLvl>& vtx_in_cell_lvl,
                             const char maskval);
  
  void get_vtx_v(const Ar2& evnmls, const int vi, const int li, Real vel[3]) const;
  Real get_vtx_z(const Ar2& ccz, const int vi, const int li) const;
  
  template <int dim>
  void interp_vel_from_vtxs_to1(const int ci, const int li, const Ar3& vtxvels,
                                const Real* p, Real* p_vel) const;
  void interp_vel_from_vtxs_to(const CellLvl& cl, const Ar3& vtxvels,
                               const int nilvl, const Real* zint, const Real* p,
                               Real* p_vel) const;
  template <int dim>
  void grad_vel_times_vel1(const int ci, const int li, const Ar3& vtxvels,
                           const Real* vel, Real* y);
  void grad_vel_times_vel(const CellLvl& cl, const Ar3& vtxvels,
                          const int nilvl, const Real* zint,
                          const Real* pos, const Real* vel, Real* y);
  template <int dim>
  void grad_vel_times_vel_at_vtx1(const int vi, const int li, const Ar3& vtxvels,
                                  const Real* vel, Real* y);
  void grad_vel_times_vel_at_vtx(const int vi, const int li, const Ar3& vtxvels,
                                 const int nilvl, const Real* vel, Real* y);
  
  int get_containing_cell_given_relations(const int* relation, const int nrelation,
                                          const Real* p);
  int get_containing_cell_given_vtx(const int vi, const Real* p);
  int get_containing_cell_given_cell(const int ci, const Real* p);
  int get_containing_level_given_z(const int nilvl, const Real* const zint,
                                   const Real z, int lvl_hint = -1);
  CellLvl get_containing_cell_lvl_given_vtx_lvl(const Levels& lmgr,
                                                const int vi, const int li,
                                                const Real* p);
  CellLvl get_containing_cell_lvl_given_cell_lvl(const Levels& lmgr,
                                                 const CellLvl& cl, const Real* p);
  int find_cell_given_starting_cell(const int ci0, const Real* const p,
                                    // if non-null, set to number of iterations
                                    int* cnt = nullptr);
  int escape_cove(const int ci, Real min_dist2, const Real* const p);

  static void
  fill_evnmls(const MpasMesh& mm, const Real time, const VelocityField& f,
              const Ar2& ecs, const Acr1& angle_edge, const Ar2& evnmls,
              const int nedge, const int nlvl, Real* const wrk = nullptr);
  static void
  fill_ccvverts(const MpasMesh& mm, const Real time, const VelocityField& f,
                const Ar2& ccs, const Ar2& ccz, const int ncell, const int nlvl);

  // For testing.
  struct TestZVelocityField;
};

#endif
