#ifndef INCLUDE_SLMMIR_REMAP_DATA_HPP
#define INCLUDE_SLMMIR_REMAP_DATA_HPP

#include "slmm_nla.hpp"

#include "slmmir.hpp"
#include "slmmir_mesh.hpp"

//   fwd = forward: The mesh at t_{n-1} is the departure mesh and is integrated
// forward in time. It is the source mesh.
//   bwd = backward: The mesh at t_n is the departure mesh and is integrated
// backward in time. It is the target mesh.
//   R = M \ T. M is the mass matrix. T is the mixed mass matrix mapping source
// to target.
class FullMassMatrix {
public:
  typedef BlockMatrix<Real, int, int> MT;

private:
  int np_;
  MT m_;

public:
  typedef std::shared_ptr<FullMassMatrix> Ptr;
  
  FullMassMatrix () : np_(0) {}
  FullMassMatrix (const int nelem, const int np) { init(nelem, np); }

  void init(const int nelem, const int np);

  int np  () const { return np_; }
  int np2 () const { return np_*np_; }
  int np4 () const { return np2()*np2(); }

  const Real* block(const int i) const;
  Real* block(const int i);

  const MT& get_M () const { return m_; }

  // Form L'L = M for all elements.
  void factor();
  // Form L'L = M.
  void factor(const int elem);
  // Solve L'L bx = bx.
  void solve(const int elem, Real* const bx, const int nrhs,
             const int ldbx) const;
  // Solve L bx = bx.
  void fwdsub(const int elem, Real* const bx, const int nrhs,
              const int ldbx) const;
  // Solve L' bx = bx.
  void bwdsub(const int elem, Real* const bx, const int nrhs,
              const int ldbx) const;

  // Solve
  //     min_x norm(A x - b) == min_x 1/2 x'A'A x - x'A'b
  //      st   c' x = d,
  // where M = A'A is the mass matrix (for one target cell), and T y = A'b is
  // the mixed mass matrix's action on source element values y (apply_T; here, T
  // is for one target cell but of course multiple source cells, so y represents
  // all possible sources).
  //   c is either c as above, or inv(L) c. If it is c, c is overwritten with
  // inv(L) c. x is T y on input and x on output.
  void solve_1eq_ls(const int elem, Real* x, const int nrhs, const int ldx,
                    Real* c, const bool c_is_Linvc, const Real d) const;
};

class RemapData {
public:
  typedef std::shared_ptr<RemapData> Ptr;
  typedef BlockMatrix<Real, Int, Int> MT;
  typedef Array<Real> VT;
  typedef siqk::Octree<SphereGeo, 10> Octree;
  struct CscGraph { Array<int> colptr, rowidx; };

  // Full block-diag target-target mass matrix, factored.
  FullMassMatrix fmm_;
  // Search tree over Eulerian mesh.
  Octree ot_;
  // Target-source matrix.
  MT T_;
  CscGraph Tgt_;
  // Jacobian(ref square -> sphere) on Eulerian mesh.
  ARealArray Jt_;
  // Eulerian mesh basis function integrals on the sphere.
  ARealArray dgbfi_, cgbfi_;

  Method::Enum method_;
  Dmc::Enum dmc_;
  // Eulerian mesh basis function integrals according to Homme definition.
  ARealArray dgbfi_gll_;
  // Pointer to either dgbfi_ or dgbfi_gll_, giving the definition of mass.
  ARealArray dgbfi_mass_;
  // Values of a block CRS matrix with same graph as T_. Each block p_ij is
  // 1xn. p_ij(k) is the proportion of source basis function k's integral in
  // overlap cell (i,j).
  Array<Real> p_s_ol_;
  // With L = chol(fmm.M) for cell e, store L \ dgbfi_mass_(idxs(e)), idxs(e)
  // the DOF indices.
  ARealArray Linv_dgbfi_mass_;

  void init_dgbfi_mass(const ARealArray& F);

public:
  RemapData(const Mesh& m, const Method::Enum method, const Dmc::Enum dmc);

  // Set up.
  MT& T () { return T_; }
  Array<Real>& p_s_ol () { return p_s_ol_; }
  CscGraph& Tgt () { return Tgt_; }

  // Apply.
  const FullMassMatrix& fmm () const { return fmm_; }
  Int T_nrows () const { return T_.M()*T_.m(); }
  Int T_ncols () const { return T_.N()*T_.n(); }
  const CscGraph& Tgt () const { return Tgt_; }
  const Octree& octree () const { return ot_; }
  const ARealArray& Jt () const { return Jt_; }
  const ARealArray& dgbfi () const { return dgbfi_; }
  const ARealArray& dgbfi_gll () const { return dgbfi_gll_; }
  const ARealArray& cgbfi () const { return cgbfi_; }
  Method::Enum method () const { return method_; }
  Dmc::Enum dmc () const { return dmc_; }
  const ARealArray& dgbfi_mass () const { return dgbfi_mass_; }
  const ARealArray& Linv_dgbfi_mass () const {
    return Linv_dgbfi_mass_;
  }

  // x = M_full \ b.
  void solve_M_full(Real* bx, const int nrhs, const int ldxb) const;

  // For eq-constrained dmc_, set up and call solve_1eq_ls.
  void remap(const Real* x, const int ldx, Real* y, const int ldy,
             const int nrhs, const Real* FsmoFtm = nullptr) const;
  // Apply remap to just cell_idx. src and tgt are the full multivectors (not
  // offset to cell ti's data).
  void remap_cell(const int cell_idx, const Real* src, const int lds,
                  Real* tgt, const int ldt, const int nrhs,
                  const Real* FsmoFtm = nullptr) const;

  // Perform and print some checks. Each entry of these Jacobians is the
  // integral over the spherical quad of a basis function. So it's really more
  // than just a Jacobian.
  void check(const Real* Js, const Real* Jt) const;
  // If T is expected to be identical to M (analytically), check how close it
  // really is. Works only before 'factor' is called.
  void compare_MT() const;
  // For analysis.
  const MT& T () const { return T_; }

private:
  inline void apply_T_cell(const int cell_idx, const Real* x, const int ldx,
                           Real* y, const int ldy, const int nrhs,
                           const Real* FsmoFtm) const;
};

#endif
