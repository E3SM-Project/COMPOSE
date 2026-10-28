#ifndef INCLUDE_MPAS_ADVECT_IR_TEST_HPP
#define INCLUDE_MPAS_ADVECT_IR_TEST_HPP

#include "mpas_advect_ir.hpp"
#include "mpas_gallery.hpp"
#include "mpas_netcdf.hpp"

#include <string>

namespace ir {

struct IcsType {
  enum Enum { one_2d, all_2d, dcmip12t11 };

  static int get_ntracer (const Enum it) {
    if (it == one_2d) return 1;
    if (it == all_2d) return 4;
    return 5;
  }
};

struct TestProblemOptions {
  static constexpr int max_n_ics = 5;
  static const InitialCondition::Shape all_ics[max_n_ics];
  
  InitialCondition::Shape ic;
  std::string flow_ic_str;
  int nstep, nsubstep, traj_nhalo;
  bool movie;
  IcsType::Enum ics_type;

  TestProblemOptions ()
    : flow_ic_str(""), nstep(0), nsubstep(-1), traj_nhalo(-1), movie(false),
      ics_type(IcsType::one_2d)
  {}
};

void set(const MpasMesh& mm, Comm& he, TestProblemOptions& o, Params& p);

struct OutputMgr {
  typedef std::shared_ptr<const OutputMgr> ConstPtr;
  NetcdfWriter::Ptr ncw;
  std::vector<NetcdfWriter::Desc> q_descs, rho_descs;
  int time_step_write_stride;
  OutputMgr () : time_step_write_stride(1) {}
};

void run_test_problem(
  const RemapFormOp::Ptr& op, TestProblemOptions& o,
  const OutputMgr::ConstPtr nw = nullptr,
  // If ! o.movie and ic, mid, or fin is provided, then fill it with the first
  // tracer's values. If o.movie, then in addition ic must be provided and have
  // the capacity to hold nstep+1 records of the first tracer.
  Real* ic = nullptr, Real* mid = nullptr, Real* fin = nullptr);

bool unittest(const MpasMesh::ConstPtr& mm, const siqp::sh::Mesh::ConstPtr& m,
              const Comm::Ptr& he, Real* data = nullptr);

void advect_mesh(
  const MpasMesh::ConstPtr& mm, const siqp::sh::Mesh::ConstPtr& m,
  TestProblemOptions& o,
  // Step for which to save output, starting with 1 and ending with o.nstep.
  const int stepno,
  // Oriented vertices of each cell, with max-ne the number of vertices.
  std::vector<Real>& xyzs_orig, std::vector<Real>& xyzs_advected);

} // namespace ir

#endif
