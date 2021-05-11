#ifndef INCLUDE_SLMMIR_PHYSGRID
#define INCLUDE_SLMMIR_PHYSGRID

#include "slmm_nla.hpp"
#include "slmm_gallery.hpp"
#include "slmm_vis.hpp"
using namespace slmm;

#include "slmmir.hpp"
#include "slmmir_d2c.hpp"
#include "slmmir_mesh.hpp"

namespace pg {

class Remap {
protected:
  Int np, nphys;
  Basis::ConstPtr b;

  // Diag of mass matrix for np: just the product of GLL weights.
  Array<Real> w_dd;

  void init(Int np, Int nphys, const Basis::ConstPtr& basis);
  static void init_M_dp(const Int np, const Int nphys, const Basis& b,
                        Array<Real>& M_dd);

public:
  Int get_np () const { return np; }
  Int get_nphys () const { return nphys; }
  void print() const;
};

class Gll2Fv : public Remap {
  Array<Real> M_pp, M_dp;

  void init(Int np, Int nphys, const Basis::ConstPtr& basis);
  void init_matrices();

public:
  typedef std::shared_ptr<Gll2Fv> Ptr;

  Gll2Fv (Int np, Int nphys, const Basis::ConstPtr& b) { init(np, nphys, b); }

  void print() const;

  // Remap a general density field.
  void remapd(const Real* gll_metdet, const Real* fv_metdet,
              const Real* d, Real* p) const;
  // Remap total mass density rho and tracer q.
  void remap(const Real* gll_metdet, const Real* fv_metdet,
             const Limiter::Enum limiter,
             const Real* drho, const Real* dq,
             Real* prho, Real* pq, const bool remap_rho=true) const;
};

struct Fv2Gll : public Remap {
  struct Type {
    // For smooth problems, choose idem; for discontinuous, elrecon. Each has at
    // least OOA 2 on C^1 problems. l2 is a reference and l2ep is an experiment.
    enum Enum {
      idem,   // OOA nphys but noisy on discontinuities
      l2,     // OOA 1 but smoothes discontinuities nicely
      l2ep,   // OOA 2 for 2 <= nphys <= 3, OOA 1 otherwise; equiv to idem for nphys = 2; smoothes
      elrecon // OOA 2 for nphys >= 2; equiv to idem for nphys = 2; a little less smoothing than l2.
    };
    static Enum convert(const std::string& s);
    static std::string convert(const Enum e);
  };

  typedef std::shared_ptr<Fv2Gll> Ptr;
  typedef std::shared_ptr<const Fv2Gll> ConstPtr;

  virtual void print() const {}
  virtual void remapd(const Real* gll_metdet, const Real* fv_metdet,
                      const Real* p, Real* d) const = 0;
  // Remap total mass density rho and tracer q.
  virtual void remap(const Real* gll_metdet, const Real* fv_metdet,
                     const Limiter::Enum limiter, const Real pqlo, const Real pqhi,
                     const Real* prho, const Real* pq,
                     Real* drho, Real* dq, const bool remap_rho=true) const;
  virtual void reconstruct_nphys1(const Real* metdet, Real* rho, Real* q) const {}

  static Fv2Gll::Ptr create(const Fv2Gll::Type::Enum type, Int np, Int nphys,
                            const Basis::ConstPtr& b);
};

class IdemFv2Gll : public Fv2Gll {
  static const Int np_nphys1;
  Int npi; // internal np
  Array<Real> op_p_to_d, npi_to_np;

  void init(Int np, Int nphys, const Basis::ConstPtr& basis);
  void init_matrices();
  void remapd(const Real* gll_metdet, const Real* fv_metdet,
              const Real* p, Real* d) const override;
  // Interp matrices are row-major.
  static void build_npi_to_np_matrix(const Int np_from, const Int np_to,
                                     const Basis& b, Array<Real>& interp);
  void reconstructd_nphys1(Real* rho) const;

public:
  IdemFv2Gll (Int np, Int nphys, const Basis::ConstPtr& b) { init(np, nphys, b); }
  void reconstruct_nphys1(const Real* metdet, Real* rho, Real* q) const override;
};

class L2Fv2Gll : public Fv2Gll {
  Int npi; // internal np
  Array<Real> op_p_to_d;

  void init(Int np, Int nphys, const Basis::ConstPtr& basis);
  void init_matrices();
  void remapd(const Real* gll_metdet, const Real* fv_metdet,
              const Real* p, Real* d) const override;

public:
  L2Fv2Gll (Int np, Int nphys, const Basis::ConstPtr& b) { init(np, nphys, b); }
};

class L2ExceptPerimFv2Gll : public Fv2Gll {
  Int npi; // internal np
  Array<Real> op_p_to_d;

  void init(Int np, Int nphys, const Basis::ConstPtr& basis);
  void init_matrices();
  void remapd(const Real* gll_metdet, const Real* fv_metdet,
              const Real* p, Real* d) const override;

public:
  L2ExceptPerimFv2Gll (Int np, Int nphys, const Basis::ConstPtr& b) { init(np, nphys, b); }
};

class ElemLclReconFv2Gll : public Fv2Gll {
  Array<Real> op_p_to_d;

  void init(Int np, Int nphys, const Basis::ConstPtr& basis);
  void init_matrices();
  void remapd(const Real* gll_metdet, const Real* fv_metdet,
              const Real* p, Real* d) const override;

public:
  ElemLclReconFv2Gll (Int np, Int nphys, const Basis::ConstPtr& b) { init(np, nphys, b); }
};

struct MeshData {
  typedef std::shared_ptr<MeshData> Ptr;

  Int nelem, ncgll, ndgll;

  // Jacobian of ref elem point -> sphere, ordered as DGLL mesh points.
  ARealArray gll_metdet;
  // metdet w(i) w(j), ordered as DGLL mesh points.
  ARealArray spheremp;
  D2Cer::Ptr d2cer;

  Array<Real> fv_metdet, wrk;

  AIdxArray geo_c2cnbrs_ptr, geo_c2cnbrs;

  MeshData (const Mesh& m, const Real* jacobian_gll) { init(m, jacobian_gll); }

private:
  void init(const Mesh& m, const Real* jacobian_gll);
};

struct PhysgridOps {
  typedef std::shared_ptr<PhysgridOps> Ptr;

  const Int np, nphys;
  Basis::ConstPtr basis;
  MeshData::Ptr mesh_data;
  Gll2Fv::Ptr gll2fv;
  Fv2Gll::Ptr fv2gll;

  PhysgridOps(const Mesh& m, const Int np, const Int nphys,
              const Fv2Gll::Type::Enum fv2gll_type = Fv2Gll::Type::idem,
              const Real* jacobian_gll = nullptr);
};

struct PhysgridData {
  typedef std::shared_ptr<PhysgridData> Ptr;

  PhysgridOps::Ptr ops;
  Limiter::Enum limiter;
  Array<Real> lat, lon, rho, q;
};

vis::VisWriter::Ptr make_pg_vis(const std::string& fname_prefix,
                                Int ne, Int nphys,
                                Int res=256, bool nonuniform=false);

} // namespace pg

#endif
