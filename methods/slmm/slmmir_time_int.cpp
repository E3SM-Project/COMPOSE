#include "slmmir_time_int.hpp"
#include "slmm_time_int.hpp"
#include "slmm_gll.hpp"
using SphereGeo = siqk::SphereGeometry;

MeshIntegrator::MeshIntegrator () : best_possible_accuracy(false) {}
void MeshIntegrator::config_for_best_possible_accuracy () { best_possible_accuracy = true; }

LatLonMeshIntegrator::LatLonMeshIntegrator (const AVec3s& p) {
  const Int nn = nslices(p);
  ll_.resize(2*nn);
# pragma omp parallel for
  for (Int i = 0; i < nn; ++i) {
    const auto n = slice(p, i);
    Real* const lli = ll_.data() + 2*i;
    xyz2ll(n[0], n[1], n[2], lli[0], lli[1]);
  }
}

template<typename OdeFn>
void MeshIntegratorWithOdeFn<OdeFn>
::integrate (const Real ts, const Real tf, AVec3s& p) {
  const Int nn = nslices(p);
  assert(2*nn <= static_cast<Int>(ll_.size()));
  ws_.resize(omp_get_max_threads());
# pragma omp parallel for schedule(static, 4)
  for (Int i = 0; i < nn; ++i) {
    const int tid = omp_get_thread_num();

    // Our primary interest in these numerical experiments is order of
    // accuracy when the flow field is exact. Hence here we use extremely
    // tight error tolerances.
    timeint::Options opts;
    opts.set_abs_tol(std::numeric_limits<Real>::epsilon());
    opts.set_rel_tol(
#ifdef RELAX_TIME
      (best_possible_accuracy ? 1e2 : 1e8)
#else
      1e2
#endif
      *std::numeric_limits<Real>::epsilon());
    opts.set_initial_step(initial_step_[i]);

    timeint::Info info;
    OdeFn fun;
    fun.set_xyz_form(use_xyz_form_);
    if ( ! use_xyz_form_) {
      Real lli[] = {ll_[2*i], ll_[2*i+1]};
      timeint::ark45(opts, fun, lli, 2, ts, tf, ws_[tid], &info);
      auto n = slice(p, i);
      ll2xyz(lli[0], lli[1], n[0], n[1], n[2]);
    } else {
      Real u[3];
      ll2xyz(ll_[2*i], ll_[2*i+1], u[0], u[1], u[2]);
      timeint::ark45(opts, fun, u, 3, ts, tf, ws_[tid], &info);
      SphereGeo::normalize(u);
      auto n = slice(p, i);
      for (Int j = 0; j < 3; ++j) n[j] = u[j];
    }
    initial_step_[i] = info.good_initial_step;
  }
}

MeshRotator::MeshRotator(const AVec3s& p)
  : LatLonMeshIntegrator(p)
{
  axis_[0] = 0.2; axis_[1] = 0.7; axis_[2] = 1;
  SphereGeo::normalize(axis_);
  resize(p_, nslices(p));
  copy(p_, p);
}

void MeshRotator
::integrate (const Real ts, const Real tf, AVec3s& p) {
  const Int nn = nslices(p);
  assert(2*nn <= static_cast<Int>(ll_.size()));
  const Real
    period = day2sec(12),
    a = 2*M_PI*(tf - ts)/period;
  Real r[9];
  form_rotation(axis_, a, r);
# pragma omp parallel for
  for (Int i = 0; i < nn; ++i) {
    auto n = slice(p_, i);
    const Real x = n[0], y = n[1], z = n[2];
    n = slice(p, i);
    n[0] = r[0]*x + r[1]*y + r[2]*z;
    n[1] = r[3]*x + r[4]*y + r[5]*z;
    n[2] = r[6]*x + r[7]*y + r[8]*z;
  }
}

std::shared_ptr<MeshIntegrator> MeshIntegratorFactory
::create (const std::string& ode, const bool use_xyz_form,
          const AVec3s& p)
{ return create(from_string(ode), use_xyz_form, p); }

std::shared_ptr<MeshIntegrator> MeshIntegratorFactory
::create (const Enum& ode, const bool use_xyz_form,
          const AVec3s& p) {
  switch (ode) {
  case Dcmip1d3ll:
    return std::make_shared<MeshIntegratorWithOdeFn<
      gallery::Dcmip1d3llOdeFn> >(p, use_xyz_form);
  case NonDivergentWindField:
    return std::make_shared<MeshIntegratorWithOdeFn<
      gallery::NonDivergentWindField> >(p, use_xyz_form);
  case DivergentWindField:
    return std::make_shared<MeshIntegratorWithOdeFn<
      gallery::DivergentWindField> >(p, use_xyz_form);
  case NonDivergentWindFieldHack:
    return std::make_shared<MeshIntegratorWithOdeFn<
      gallery::NonDivergentWindFieldHack> >(p, use_xyz_form);
  case Rotate:
    return std::make_shared<MeshRotator>(p);
  case MovingVortices:
    return std::make_shared<MeshIntegratorWithOdeFn<
      gallery::MovingVortices> >(p, use_xyz_form);
  case TestWindField:
    return std::make_shared<MeshIntegratorWithOdeFn<
      gallery::TestWindField> >(p, use_xyz_form);
  default:
    assert(0);
    return nullptr;
  }
}

VelocityInterpolatorMeshIntegrator::VelocityInterpolatorMeshIntegrator (
  const MeshIntegrator::Ptr& imi_coarse, const Mesh::Ptr& im_coarse,
  const Mesh::Ptr& im_fine)
  : mi_coarse(imi_coarse), m_coarse(im_coarse), m_fine(im_fine), p_coarse(m_coarse->cgll_p),
    // Move the mesh using GLL basis.
    mesh_interp(m_coarse->np, Basis::create(Basis::Type::gll), m_fine->np, m_fine->basis)
{
  // The coarse mesh should always have GLL nodes.
  assert(m_coarse->basis->gll_nodes());
  ne = nslices(m_coarse->dglln2cglln) / square(m_coarse->np);
  assert(ne == nslices(m_fine->dglln2cglln) / square(m_fine->np));
}

void VelocityInterpolatorMeshIntegrator::
integrate (const Real ts, const Real tf, AVec3s& p_fine) {
  assert(nslices(p_fine) == nslices(m_fine->cgll_p));
  mi_coarse->integrate(ts, tf, p_coarse);
  mesh_interp.apply(*m_coarse, *m_fine, p_coarse, p_fine);
}
