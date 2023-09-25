#ifndef INCLUDE_SLMMIR_REMAPPER_HPP
#define INCLUDE_SLMMIR_REMAPPER_HPP

#include "slmm_fit_extremum.hpp"

#include "slmmir_p_refine.hpp"
#include "slmmir_physgrid.hpp"
#include "slmmir_mono_data.hpp"
#include "slmmir_d2c.hpp"
#include "slmmir.hpp"

struct C2DRelations { AIdxArray ptr, r; };

// Orchestrate the remap operation in one time step.
//   "np basis" is the full nonmonotonic basis.
//   "npmono basis" is the monotonic basis used in the cell mean boundedness
// correction.
//   cmbc: cell mean boundedness correction.
class Remapper {
  const Int np_, np2_, ncell_, dnn_, cnn_;

  const std::shared_ptr<const Mesh> m_;
  const std::shared_ptr<RemapData> rd_;
  const std::shared_ptr<const MonoData> md_;
  spf::MassRedistributor::Ptr mrd_;
  D2Cer::Ptr d2cer_;

  // Basis function integral on the source density field's mesh, divided by the
  // integral on the target mesh. In IR (as opposed to CDG), this quotient is
  // the density factor. Intuitively, the idea is to multiply a nodal density
  // value by its basis function to get a mass, then divide it by the basis
  // function on the target mesh to get the new density.
  ARealArray FsmoFtm_;

  // Density in various states of remap. rho_src_* do not have the the Jt/Js
  // density update because they live on the Eulerian mesh.
  ARealArray
    rho_src_cmbc_,      // Corrected, in the np basis.
    rho_tgt_cmbc_;      // Density in target cell given by cmbc sources.

  // For ISL.
public:
  class IslImpl;
private:
  std::shared_ptr<IslImpl> isl_impl_;
  C2DRelations c2d_v_, c2d_f_;
  AIdxArray cn_src_cell_;
  ARealArray2 q_data_;
  ARealArray Je_, Jdep_;

  Real perturb_rho_;
  bool record_total_mass_redistribution_;
  std::vector<Real> total_mass_redistribution_, total_mass_discrepancy_;

  bool subcell_bounds_;

  std::shared_ptr<FitExtremum> fit_extremum_;

  PRefineData::Ptr pr_;

private:
  void dgll_set_src_tgt_in(Real* const src_in, Real* const tgt_in,
                           Real*& src, Real*& tgt);

  void dgll_set_tgt_out(const Real* const tgt, Real* const tgt_out);

  void project_nolimiter(
    Real* const src_tracer, Real* const tgt_tracer, const Int ntracers,
    Real* const src_density, Real* const tgt_density);

  void apply_R(const Real* const src, Real* const tgt) const;

  void apply_R(const Real* const src, Real* const tgt,
               const Int* const tis, const Int nti) const;

  void perturb_rho(Real* const rho, const Real p);

  void project_and_limit_cdr(
    Real* const src_tracer, Real* const tgt_tracer, const Int ntracers,
    Real* const src_density, Real* const tgt_density,
    const Filter::Enum cdr_method, const Real rho_perturbation);

  void project_and_limit(
    Real* const src_tracer, Real* const tgt_tracer, const Int ntracers,
    Real* const src_density, Real* const tgt_density);

  void init_isl();
  void init_isl_jacobian(const Mesh& m);
  static void make_c2d_relations(
    const Int cnn, const AIdxArray& dglln2cglln,
    C2DRelations& c2d);

  void interp(
    const Mesh& m, const C2DRelations& c2d, const AVec3s& advected_p,
    const Real* const src_tracer, Real* const tgt_tracer, const Int ntracers,
    const Real* const src_density, Real* const tgt_density,
    const bool rho_isl = false);

  // Constrained density reconstruction for classical semi-Lagrangian.
  void isl_cdr(
    const Mesh& m, const RemapData& rd_src, const RemapData& rd_tgt,
    const MonoData& md, const spf::MassRedistributor::Ptr& mrd,
    Real* const src_density, Real* const src_tracer, const Int np_src,
    Real* const tgt_density, Real* const tgt_tracer, const Int ntracers,
    const bool positive_only, const bool rho_isl, const Filter::Enum cdr_method);
  void isl_cdr_rho(
    const Mesh& m, const RemapData& rd_src, const RemapData& rd_tgt,
    const MonoData& md, const spf::MassRedistributor::Ptr& mrd,
    Real* const src_density, const Int np2_src,
    Real* const tgt_density, const Int np2_tgt,
    const Int ntracers);

public:
  Remapper(const PRefineData::Ptr& pr, const MonoData::Ptr& md,
           const D2Cer::Ptr& d2cer);

  void use_subcell_bounds();
  void use_fit_extremum();

  // Stuff to analyze tracer consistency via project_and_limit_cdr.
  void set_rho_perturbation(const Real p);
  void record_total_mass_redistribution(const bool record);
  void track_footprint(const bool track);
  Real get_total_redistributed_mass(const Int& tracer_idx);
  Real get_total_mass_discrepancy(const Int& tracer_idx);

  // On input, src_tracer is rho*tracer. On output, it is just the updated
  // tracer. Density is removed for output and error checking.
  void remap(const AVec3s& advected_p,
             Real* const src_tracer, Real* const tgt_tracer, const Int ntracers,
             Real* const src_density, Real* const tgt_density);
    
  void isl(const AVec3s& advected_p, Real* const src_tracer,
           Real* const tgt_tracer, const Int ntracers,
           Real* const src_density, Real* const tgt_density,
           const bool positive_only, const bool rho_isl,
           const Real* toychem_tendency, Int toychem_idx_beg, Int toychem_idx_end,
           const pg::PhysgridData::Ptr& pg);

  void init_physgrid(const Real* const rho, const Real* const tracer,
                     const Int ntracers, pg::PhysgridData& pg);

  // For p-refined ISL analysis only.
  void get_prefine_internal(const Real*& q, const Real*& dgbfi_mass, Int& len,
                            const bool first) const;
};

#endif
