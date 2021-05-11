#ifndef INCLUDE_SLMMIR_TIME_INT_HPP
#define INCLUDE_SLMMIR_TIME_INT_HPP

#include "slmm_time_int.hpp"
#include "slmm_gallery.hpp"
using namespace slmm;

#include "slmmir_mesh.hpp"
#include "slmmir_p_refine.hpp"

#include <memory>
#include <vector>

struct TimeInt {
  enum Enum { exact, line, interp, interpline };
  static std::string convert (const Enum& e) {
    switch (e) {
    case Enum::exact: return "exact";
    case Enum::line: return "line";
    case Enum::interp: return "interp";
    case Enum::interpline: return "interpline";
    default: throw std::runtime_error("TimeInt: Not a valid enum.");
    }
  }
  static Enum convert (const std::string& s) {
    if (eq(s, "exact")) return Enum::exact;
    if (eq(s, "line")) return Enum::line;
    if (eq(s, "interp")) return Enum::interp;
    if (eq(s, "interpline")) return Enum::interpline;
    throw std::runtime_error(std::string("TimeInt: Not a valid string: ") + s);
  }
  static bool is_interp (const Enum& e) { return e == interp || e == interpline; }
};

struct MeshIntegrator {
  typedef std::shared_ptr<MeshIntegrator> Ptr;
  MeshIntegrator();
  virtual ~MeshIntegrator () {}
  virtual Int nnodes() const = 0;
  virtual void integrate(const Real ts, const Real tf, AVec3s& p) =0;
  // For almost all tests of interest, this is not needed and causes substantial
  // performance loss. But for the midpoint check or any other very long
  // trajectory integration, we want to boost accuracy as much as possible. Call
  // this to do that.
  void config_for_best_possible_accuracy();
protected:
  bool best_possible_accuracy;
};

class LatLonMeshIntegrator : public MeshIntegrator {
protected:
  std::vector<Real> ll_;
public:
  typedef std::shared_ptr<MeshIntegrator> Ptr;
  LatLonMeshIntegrator(const AVec3s& p);
  virtual ~LatLonMeshIntegrator () {}
  virtual Int nnodes () const override { return ll_.size() / 2; }
  virtual void integrate(const Real ts, const Real tf, AVec3s& p) =0;
};

template<typename OdeFn>
class MeshIntegratorWithOdeFn : public LatLonMeshIntegrator {
  std::vector<timeint::Workspace> ws_;
  std::vector<Real> initial_step_;
  bool use_xyz_form_;

public:
  MeshIntegratorWithOdeFn (const AVec3s& p,
                           const bool use_xyz_form = false)
    : LatLonMeshIntegrator(p), initial_step_(nslices(p), 1e-3),
      use_xyz_form_(use_xyz_form)
  {}

  void integrate(const Real ts, const Real tf, AVec3s& p) override;
};

class MeshRotator : public LatLonMeshIntegrator {
  AVec3s p_;
  Real axis_[3];

public:
  MeshRotator(const AVec3s& p);

  void integrate(const Real ts, const Real tf, AVec3s& p) override;
};

struct MeshIntegratorFactory : public gallery::WindFieldType {
  static std::shared_ptr<MeshIntegrator>
  create(const std::string& ode, const bool use_xyz_form,
         const AVec3s& p);

  static std::shared_ptr<MeshIntegrator>
  create(const Enum& ode, const bool use_xyz_form,
         const AVec3s& p);
};

struct StudyTimeIntegratorOptions {
};

struct StudyTimeIntegratorFactory : public gallery::WindFieldType {
  static std::shared_ptr<MeshIntegrator>
  create(const std::string& ode, const StudyTimeIntegratorOptions& o,
         const AVec3s& p);

  static std::shared_ptr<MeshIntegrator>
  create(const Enum& ode, const StudyTimeIntegratorOptions& o,
         const AVec3s& p);
};

struct VelocityInterpolatorMeshIntegrator : public MeshIntegrator {
  VelocityInterpolatorMeshIntegrator(const MeshIntegrator::Ptr& mi_coarse,
                                     const Mesh::Ptr& m_coarse, const Mesh::Ptr& m_fine);
  Int nnodes() const override { return nslices(m_fine->cgll_p); }
  void integrate(const Real ts, const Real tf, AVec3s& p) override;
private:
  MeshIntegrator::Ptr mi_coarse;
  Mesh::Ptr m_coarse, m_fine;
  AVec3s p_coarse;
  MeshInterpolator mesh_interp;
  Int ne;
};

#endif
