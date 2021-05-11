#include "slmmir_time_int.hpp"
#include "slmm_time_int.hpp"
using SphereGeo = siqk::SphereGeometry;

std::shared_ptr<gallery::OdeFn>
create (const gallery::WindFieldType::Enum& ode) {
  switch (ode) {
  case gallery::WindFieldType::Dcmip1d3ll:
    return std::make_shared<gallery::Dcmip1d3llOdeFn>();
  case gallery::WindFieldType::NonDivergentWindField:
    return std::make_shared<gallery::NonDivergentWindField>();
  case gallery::WindFieldType::DivergentWindField:
    return std::make_shared<gallery::DivergentWindField>();
  case gallery::WindFieldType::NonDivergentWindFieldHack:
    return std::make_shared<gallery::NonDivergentWindFieldHack>();
  case gallery::WindFieldType::MovingVortices:
    return std::make_shared<gallery::MovingVortices>();
  case gallery::WindFieldType::Rotate:
    throw std::runtime_error("Rotation field is not currently supported.");
  default:
    throw std::runtime_error("Not a WindFieldType.");
    return nullptr;
  }
}

class StudyTimeIntegrator : public LatLonMeshIntegrator {
  std::shared_ptr<gallery::OdeFn> odefn_;

public:
  StudyTimeIntegrator (
    const AVec3s& p, const gallery::WindFieldType::Enum& ode)
    : LatLonMeshIntegrator(p)
  {
    odefn_ = create(ode);
    odefn_->set_xyz_form(true);
  }

  void integrate(const Real ts, const Real tf, AVec3s& p) override {
    const Int nn = nslices(p);
    assert(2*nn == static_cast<Int>(ll_.size()));
    const Real dt = tf - ts;
#   pragma omp parallel for
    for (Int i = 0; i < nn; ++i) {
      Real u[3];
      ll2xyz(ll_[2*i], ll_[2*i+1], u[0], u[1], u[2]);

      Real uh[3], f[3];
      for (int j = 0; j < 3; ++j) uh[j] = u[j];
      const Real th = 0.5*(ts + tf);
      for (int it = 0; it < 2; ++it) {
        odefn_->eval(th, uh, f);
        for (int j = 0; j < 3; ++j) uh[j] = u[j] + 0.5*dt*f[j];
      }
      for (int j = 0; j < 3; ++j) u[j] += dt*f[j];
      
      // Normalization doesn't matter. SphereGeo::normalize(u);
      auto n = slice(p, i);
      for (Int j = 0; j < 3; ++j) n[j] = u[j];
    }
  }
};

std::shared_ptr<MeshIntegrator> StudyTimeIntegratorFactory
::create (const std::string& ode, const StudyTimeIntegratorOptions& o,
          const AVec3s& p) {
  return create(from_string(ode), o, p);
}

std::shared_ptr<MeshIntegrator> StudyTimeIntegratorFactory
::create (const gallery::WindFieldType::Enum& ode,
          const StudyTimeIntegratorOptions& o,
          const AVec3s& p) {
  return std::make_shared<StudyTimeIntegrator>(p, ode);
}
