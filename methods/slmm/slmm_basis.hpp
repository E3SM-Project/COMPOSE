#ifndef INCLUDE_SLMM_BASIS_HPP
#define INCLUDE_SLMM_BASIS_HPP

#include "slmm_defs.hpp"

#include <memory>

namespace slmm {

struct Basis {
  typedef std::shared_ptr<Basis> Ptr;
  typedef std::shared_ptr<const Basis> ConstPtr;

  enum { np_max = SLMM_NP_MAX };

  struct Type {
    enum Enum { gll, gll_offset_nodal, gll_nodal, free_nodal, constant_cell,
                uniform_offset_nodal, uniform_reduced };
    static Enum convert(const std::string& s);
    static std::string convert(const Enum e);
  };

  static Basis::Ptr create(const Type::Enum e);
  static Basis::Ptr create_basis_from_string(const std::string& s);

  // Return false if there is a failure, such as unsupported np.
  virtual ~Basis () {}
  virtual bool gll_nodes () const { return true; }
  virtual const char* name () const { return "Basis"; }
  virtual bool get_x(const Int& np, const Real*& coord) const = 0;
  virtual bool get_w(const Int& np, const Real*& wt) const = 0;
  // Max polynomial degree in basis.
  virtual Int max_degree(const Int& np) const { return np-1; }
  virtual bool eval(const Int& np, const Real& x, Real* const v) const = 0;
  virtual bool eval_derivative(const Int& np, const Real& x, Real* const v) const = 0;

  template <typename Scalar>
  static KOKKOS_INLINE_FUNCTION
  void eval_lagrange_poly (const Int& np, const Scalar* const x_gll,
                           const Scalar& x, Scalar* const y) {
    for (int i = 0; i < np; ++i) {
      Scalar f = 1;
      for (int j = 0; j < np; ++j)
        f *= (j == i ?
              1 :
              (x - x_gll[j]) / (x_gll[i] - x_gll[j]));
      y[i] = f;
    }
  }

  template <typename Scalar>
  static KOKKOS_INLINE_FUNCTION
  void eval_lagrange_poly_derivative (const Int& np, const Scalar* const x_gll,
                                      const Scalar& x, Scalar* const y) {
    for (int i = 0; i < np; ++i) {
      Scalar f = 0;
      for (int j = 0; j < np; ++j) {
        if (j == i) continue;
        Scalar g = 1;
        for (int k = 0; k < np; ++k)
          g *= ((k == i) ? 1 :
                ((k == j ? 1 : (x - x_gll[k])) /
                 (x_gll[i] - x_gll[k])));
        f += g;
      }
      y[i] = f;
    }
  }

  static void compute_weights(const Basis& basis, const Int np, Real* wts);
  static Int compute_and_print_weights(const Basis& basis, bool print_x = false,
                                       bool test = false);

  // For physgrid, M[i*np^2 + j] is the integral of 2D basis function j over 2D
  // subcell i, with 2D basis functions and subcells ordered as in slmm::mesh.
  static bool calc_common_refinement(const Basis& b, const Int np,
                                     const Int nphys, std::vector<Real>& xcom,
                                     const bool test = false);
  static bool compute_integrals_over_subcells_2d(const Basis& b, const Int np,
                                                 const Int nphys, Real* M,
                                                 const bool test = false);

  // The row-major matrix M has range br and domain bc.
  static bool calc_common_refinement(const Basis& br, const Int npr,
                                    const Basis& bc, const Int npc,
                                    std::vector<Real>& xcom, const bool test = false);
  static bool compute_mass_matrix_2d(const Basis& br, const Int npr,
                                     const Basis& bc, const Int npc, Real* M,
                                     const bool test = false);
};

} // namespace slmm

#endif
