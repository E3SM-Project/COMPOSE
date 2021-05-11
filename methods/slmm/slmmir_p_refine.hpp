#ifndef INCLUDE_SLMMIR_P_REFINE_HPP
#define INCLUDE_SLMMIR_P_REFINE_HPP

#include <memory>

#include "slmmir_mono_data.hpp"

struct PRefineData {
  typedef std::shared_ptr<PRefineData> Ptr;

  bool prefine;
  int experiment;
  // Velocity, fine, and transport meshes. Transport mesh could be either one,
  // depending on p_refine_experiment.
  Mesh::Ptr m_v, m_f, m_t;
  RemapData::Ptr rd_v, rd_f, rd_t;
  MonoData::Ptr md_v, md_f, md_t;
};

class Transferer2D {
protected:
  Int np_from, np_to, np_from2, np_to2;
  Basis::Ptr b_from, b_to;

public:
  typedef std::shared_ptr<Transferer2D> Ptr;

  Transferer2D(Int np_from, const Basis::Ptr& b_from, Int np_to, const Basis::Ptr& b_to);
  ~Transferer2D () {}

  virtual void apply(const Real* src, Real* dst) const = 0;

  Int get_np_from () const { return np_from; }
  Int get_np_to () const { return np_to; }
  const Basis& get_basis_from () const { return *b_from; }
  const Basis& get_basis_to () const { return *b_to; }
};

class MatrixBasedTransferer2D : public Transferer2D {
protected:
  Array<Real> op;

public:
  MatrixBasedTransferer2D(Int np_from, const Basis::Ptr& b_from, Int np_to, const Basis::Ptr& b_to);
  virtual void apply(const Real* src, Real* dst) const override;
};

struct Interpolator2D : public MatrixBasedTransferer2D {
  typedef std::shared_ptr<Interpolator2D> Ptr;

  Interpolator2D(Int np_from, const Basis::Ptr& b_from, Int np_to, const Basis::Ptr& b_to);

  /* Solve
        x = arg min_x 1/2 norm(x - y, 2)^2
                 st   B x = u,
     where B is op in this class. The Lagrangian is
        L := 1/2 norm(x - y, 2)^2 + lam'(B x - u),
     so
        x - y + B'lam = 0
              B x - u = 0.
     The expected setting is that B interpolates a higher order to a lower one:
     find the high-order field x nearest y that interpolates u.
  */
  void find_x_nearest_y_interpolating_u(const Real* y, const Real* u, Real* x);

private:
  QrFac::Ptr qr_fac;
  
  void init();
};

class MeshInterpolator {
  Interpolator2D i2d;
public:
  typedef std::shared_ptr<MeshInterpolator> Ptr;
  MeshInterpolator (Int np_from, const Basis::Ptr& b_from, Int np_to, const Basis::Ptr& b_to)
    : i2d(np_from, b_from, np_to, b_to) {}
  void apply(const Mesh& m_coarse, const Mesh& m_fine,
             const AVec3s& p_coarse,
             AVec3s& p_fine) const;
  const Interpolator2D& get_i2d () const { return i2d; }
};

void calc_pref_gll_quantities(const Mesh& mc, const RemapData& rdc, const Mesh& mf,
                              RemapData& rdf);

#endif
