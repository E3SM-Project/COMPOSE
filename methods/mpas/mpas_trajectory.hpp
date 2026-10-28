#ifndef INCLUDE_MPAS_TRAJECTORY_HPP
#define INCLUDE_MPAS_TRAJECTORY_HPP

#include "mpas_mesh.hpp"
#include "mpas_siqp.hpp"
#include "mpas_gallery.hpp"
#include "mpas_comm.hpp"
#include "mpas_time_int.hpp"

#include <memory>

// All inputs must be for unit-radius sphere, and all outputs are.
struct Trajectory {
  typedef std::shared_ptr<Trajectory> Ptr;
  typedef std::shared_ptr<const Trajectory> ConstPtr;

  struct Method {
    enum Enum { exact,      // adaptive RK45 with all vel components at vtxs
                order1, order2, taylor2 // algorithms available in practice
    };
    static Method::Enum from_short(const char sym);
    static const char* to_string(Enum e);
    static bool is_exact (const Enum e) { return e == exact; }
  };

  struct OdeFn : public timeint::OdeFn {
    const VelocityField& vf;
    const int lvl;
    const bool zcomp;
    OdeFn(const VelocityField& vf_, const int lvl_);
    bool eval(Real t, const Real* y, Real* f) const override;
  };

  struct CellLvl {
    CellLvl () : ci(-1), lvl(-1) {}
    CellLvl (const int ci_, const int lvl_) : ci(ci_), lvl(lvl_) {}
    int cell () const { return ci; }
    int level () const { return lvl; }
  private:
    int ci, lvl;
  };

  Trajectory(const MpasMesh::ConstPtr& m, const Comm::Ptr& he, VelocityField::Ptr vf,
             const Method::Enum method, const int nsubstep, const int nhalo = max_nhalo);

  void set_method_exact_rel_tol(const Real reltol);

  int num_levels() const;

  size_t get_rwrk_size() const;

  // Advect vertices backward in time from ts+dt to ts.
  void advect_vtxs(const Array2D<const Real>& p, const Real ts, const Real dt,
                   const Ar1& rwrk, const Array3D<Real>& p_advect);

  bool has_in_cell_data() const;
  // Get Eulerian cell local ID and midpoint level containing the advectect
  // vertex that started at Eulerian vertex vi and level vlvl.
  const CellLvl& get_cell_lvl_containing_vtx(const int vi, const int lvl) const;

  bool unittest(const Array2D<const Real>& vtxs, Real* data = nullptr);
  static bool convtest(const MpasMesh::ConstPtr& m, const Comm::Ptr& he,
                       Real* data = nullptr);

  static int test_search(const MpasMesh::ConstPtr& m, const Comm::Ptr& he,
                         const bool print_stats = false);

public: // for internal use
  struct Impl;
  std::shared_ptr<Impl> impl;
};

// Adjust nsubstep and nhalo based on he.szhalo() if he.parallel() and the
// method depends on halo exchanges.
void adjust_nsubstep_nhalo_if_needed(
  const Comm& he, const Trajectory::Method::Enum method,
  // Number of extra cell haloes needed at my proc's edges that cannot be source
  // cells, needed e.g. for high-order reconstruction.
  const int comm_edge_cell_halo_extra,
  int& nsubstep, int& traj_nhalo);

bool unittest_gallery(const MpasMesh::Ptr& mm, const siqp::sh::Mesh& m,
                      const Comm::Ptr& he);

#endif
