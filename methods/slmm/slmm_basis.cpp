#include "slmm_basis.hpp"
#include "slmm_gll.hpp"
#include "slmm_islet.hpp"
#include "slmm_basis_reduced.hpp"
#include "slmm_util.hpp"

#undef NDEBUG
#include <cassert>

#include <set>

namespace slmm {

Basis::Type::Enum Basis::Type::convert (const std::string& s) {
  if (s == "gll") return Type::gll;
  if (s == "gll_offset_nodal") return Type::gll_offset_nodal;
  if (s == "gll_nodal") return Type::gll_nodal;
  if (s == "free_nodal") return Type::free_nodal;
  if (s == "uniform_offset_nodal") return Type::uniform_offset_nodal;
  if (s == "uniform_reduced") return Type::uniform_reduced;
  SIQK_THROW_IF(true, "Basis::convert: not impled for " << s);
}

std::string Basis::Type::convert (const Basis::Type::Enum e) {
  switch (e) {
  case Type::gll: return "gll";
  case Type::gll_offset_nodal: return "gll_offset_nodal";
  case Type::gll_nodal: return "gll_nodal";
  case Type::free_nodal: return "free_nodal";
  case Type::uniform_offset_nodal: return "uniform_offset_nodal";
  case Type::uniform_reduced: return "uniform_reduced";
  default: SIQK_THROW_IF(true, "Basis::convert: not impled for " << e);
  }  
}

Basis::Ptr Basis::create (const Basis::Type::Enum e) {
  switch (e) {
  case Type::gll: return std::make_shared<GLL>();
  case Type::gll_offset_nodal: return std::make_shared<islet::GllOffsetNodal>();
  case Type::gll_nodal: return std::make_shared<islet::GllNodal>();
  case Type::free_nodal: return std::make_shared<islet::FreeNodal>();
  case Type::uniform_offset_nodal: return std::make_shared<islet::UniformOffsetNodal>();
  case Type::uniform_reduced: return std::make_shared<UniformNodeReduced>();
  default: SIQK_THROW_IF(true, "Basis::create: not impled for " << e);
  }
}

Basis::Ptr Basis::create_basis_from_string (const std::string& s) {
  const auto gll_nodal = s.find("x") == std::string::npos;
  if (gll_nodal) {
    const auto b = std::make_shared<islet::GllNodalFromString>();
    if ( ! b->init(s)) return nullptr;
    return b;
  } else {
    const auto b = std::make_shared<islet::FreeNodalFromString>();
    if ( ! b->init(s)) return nullptr;
    return b;
  }
}

void Basis::compute_weights (const Basis& basis, const Int np, Real* integral) {
  const Int qn = 10;
  const Real* qx, * qwt;
  GLL gll;
  gll.get_coef(qn, qx, qwt);
  const Real* xnode;
  const bool ok = basis.get_x(np, xnode);
  assert(ok);
  // Integrate numerically over each nodal region, using np=10 GLL quadrature,
  // which is enough for exact integration of degree-17 polynomials.
  Real v[np_max];
  for (Int i = 0; i < np; ++i) integral[i] = 0;
  for (Int ireg = 0; ireg < np-1; ++ireg) {
    Real reg_integral[np_max] = {0};
    for (Int qi = 0; qi < qn; ++qi) {
      const auto alpha = 0.5*(qx[qi] + 1);
      const auto x = (1 - alpha)*xnode[ireg] + alpha*xnode[ireg+1];
      basis.eval(np, x, v);
      for (Int i = 0; i < np; ++i)
        reg_integral[i] += qwt[qi]*v[i];
    }
    const auto fac = 0.5*(xnode[ireg+1] - xnode[ireg]);
    for (Int i = 0; i < np; ++i)
      integral[i] += fac*reg_integral[i];
  }
  // Numerically symmetrize.
  for (Int ireg = 0; ireg < np/2; ++ireg) {
    // std::min is to prevent a spurious -Warray-bounds warning.
    const Int other = std::min(np_max-1, np-ireg-1);
    integral[ireg] = integral[other] =
      0.5*(integral[ireg] + integral[other]);
  }
  Real total = 0;
  for (Int i = 0; i < np; ++i) total += integral[i];
  assert(reldif(2, total) < 5*std::numeric_limits<Real>::epsilon());  
  // Numerically make 2.
  const Real fac = 2/total;
  for (Int i = 0; i < np; ++i) integral[i] *= fac;
}

Int Basis::compute_and_print_weights (const Basis& basis, bool print_x, bool test) {
  const GLL* bgll = nullptr;
  Int nerr = 0;
  if (test) {
    bgll = dynamic_cast<const GLL*>(&basis);
    assert(bgll);
  }
  const Real* xnode;
  for (Int np = 2; np <= Basis::np_max; ++np) {
    // Basis need not support every np.
    if ( ! basis.get_x(np, xnode)) continue;
    Real integral[np_max] = {0};
    compute_weights(basis, np, integral);
    if ( ! test) {
      // Print C code suitable for adding to the Basis class.
      const bool odd = np % 2 == 1;
      if (print_x) {
        printf("const Real x_np%d[%d] = {-1", np, np);
        for (Int i = 1; i < np-1; ++i) {
          if (odd && i == np/2)
            printf(", 0");
          else {
            if (i < np/2) printf(", %23.16e", xnode[i]);
            else printf(", %22.16e", xnode[i]);
          }
        }
        printf(", 1};\n");
      };
      printf("const Real w_np%d[%d] = {%23.16e", np, np, integral[0]);
      for (Int i = 1; i < np; ++i) printf(", %22.16e", integral[i]);
      printf("};\n");
    }
    if (bgll) {
      // Test the numerical procedure against the known GLL weights in the case
      // that testing is requested.
      const Real* wt;
      bgll->get_w(np, wt);
      const Real tol = (np <= 13 ? 15 : 20)*std::numeric_limits<Real>::epsilon();
      for (Int i = 0; i < np; ++i)
        if (reldif(wt[i], integral[i]) > tol) {
          printf("%22.16e vs %22.16e: %9.3e\n", wt[i], integral[i],
                 reldif(wt[i], integral[i]));
          ++nerr;
        }
    }
  }
  return nerr;
}

bool Basis::calc_common_refinement (const Basis& b, const Int np,
                                    const Int nphys, std::vector<Real>& xcom,
                                    const bool test) {
  std::set<Real> s;
  const Real* x_np;
  const bool ok = b.get_x(np, x_np);
  assert(ok);
  if ( ! ok) return false;
  for (Int i = 0; i < np; ++i) s.insert(x_np[i]);
  for (Int i = 1; i < nphys; ++i) s.insert(2*Real(i)/nphys - 1);
  xcom.clear();
  for (const auto e : s) xcom.push_back(e);
  std::sort(xcom.begin(), xcom.end());
  if (test && b.gll_nodes()) {
    Int n = np + nphys - 1;
    if (np % 2 == 1 && nphys % 2 == 0) n--;
    if ((Int) xcom.size() != n) return false;
  }
  return true;
}

bool Basis
::compute_integrals_over_subcells_2d (const Basis& b, const Int np,
                                      const Int nf, Real* M, const bool test) {
  const auto subcell_idx = [&] (const Real r) -> Int {
    const Real p = (1 + r)/2;
    return std::max(0, std::min(nf-1, Int(nf*p)));
  };

  std::vector<Real> xcom;
  if ( ! calc_common_refinement(b, np, nf, xcom, test)) {
    assert(0);
    return false; 
  }
  const Int ncom = xcom.size(), np2 = square(np), nf2 = square(nf);

  Int qn = (b.max_degree(np) + 4)/2;
  if (qn > 13 && qn < 16) qn = 16;
  const Real* qx, * qwt;
  GLL gll;
  gll.get_coef(qn, qx, qwt);

  for (Int i = 0; i < np2*nf2; ++i) M[i] = 0;

  Real vx[GLL::np_max], vy[GLL::np_max];
  for (Int iy = 0; iy < ncom-1; ++iy) {
    const Real yb = xcom[iy], ye = xcom[iy+1];
    const Int ify = subcell_idx(0.5*(yb + ye));
    for (Int ix = 0; ix < ncom-1; ++ix) {
      const Real xb = xcom[ix], xe = xcom[ix+1];
      const Int ifx = subcell_idx(0.5*(xb + xe));
      Real* Msub = M + np2*(ify*nf + ifx);
      const Real fac = 0.25*(ye - yb)*(xe - xb);
      for (Int jy = 0; jy < qn; ++jy) {
        Real alpha = 0.5*(qx[jy] + 1);
        const Real y = (1 - alpha)*yb + alpha*ye;
        b.eval(np, y, vy);
        for (Int jx = 0; jx < qn; ++jx) {
          alpha = 0.5*(qx[jx] + 1);
          const Real x = (1 - alpha)*xb + alpha*xe;
          b.eval(np, x, vx);
          const Real fac1 = fac*qwt[jy]*qwt[jx];
          for (Int i = 0; i < np; ++i)
            for (Int j = 0; j < np; ++j)
              Msub[np*i + j] += fac1*vy[i]*vx[j];
        }
      }
    }
  }

  return true;
}

bool Basis::calc_common_refinement (const Basis& br, const Int npr,
                                    const Basis& bc, const Int npc,
                                    std::vector<Real>& xcom, const bool test) {
  std::set<Real> s;
  const auto insert = [&] (const Basis& b, const Int np) {
    const Real* x_np;
    const bool ok = b.get_x(np, x_np);
    assert(ok);
    for (Int i = 0; i < np; ++i) s.insert(x_np[i]);
    return true;
  };
  insert(br, npr);
  insert(bc, npc);
  xcom.clear();
  for (const auto e : s) xcom.push_back(e);
  std::sort(xcom.begin(), xcom.end());
  if (test) {
    const Int n = xcom.size();
    if (br.gll_nodes() && bc.gll_nodes()) {
      const bool e = ((npr == npc && n != npr) ||
                      (npr != npc && (npr % 2 == 1 && npc % 2 == 1) && n != npr + npc - 3) ||
                      (npr != npc && (npr % 2 == 0 || npc % 2 == 0) && n != npr + npc - 2));
      if (e) return false;
    } 
  }
  return true;
}

bool Basis::compute_mass_matrix_2d (const Basis& br, const Int npr,
                                    const Basis& bc, const Int npc, Real* M,
                                    const bool test) {
  const Real* rx, * cx;
  br.get_x(npr, rx);
  bc.get_x(npc, cx);

  std::vector<Real> xcom;
  if ( ! calc_common_refinement(br, npr, bc, npc, xcom, test)) {
    assert(0);
    return false;
  }
  const Int ncom = xcom.size(), npr2 = square(npr), npc2 = square(npc);

  Int qn = (br.max_degree(npr) + bc.max_degree(npc) + 4)/2;
  if (qn > 16) return false;
  if (qn > 13 && qn < 16) qn = 16;

  const Real* qx, * qwt;
  GLL gll;
  gll.get_coef(qn, qx, qwt);

  for (Int i = 0; i < npr2*npc2; ++i) M[i] = 0;

  Real rvx[GLL::np_max], rvy[GLL::np_max], cvx[GLL::np_max], cvy[GLL::np_max];
  for (Int iy = 0; iy < ncom-1; ++iy) {
    const Real yb = xcom[iy], ye = xcom[iy+1];
    for (Int ix = 0; ix < ncom-1; ++ix) {
      const Real xb = xcom[ix], xe = xcom[ix+1];
      const Real fac = 0.25*(ye - yb)*(xe - xb);
      for (Int jy = 0; jy < qn; ++jy) {
        Real alpha = 0.5*(qx[jy] + 1);
        const Real y = (1 - alpha)*yb + alpha*ye;
        br.eval(npr, y, rvy);
        bc.eval(npc, y, cvy);
        for (Int jx = 0; jx < qn; ++jx) {
          alpha = 0.5*(qx[jx] + 1);
          const Real x = (1 - alpha)*xb + alpha*xe;
          br.eval(npr, x, rvx);
          bc.eval(npc, x, cvx);
          const Real fac1 = fac*qwt[jy]*qwt[jx];
          for (Int ri = 0; ri < npr; ++ri)
            for (Int rj = 0; rj < npr; ++rj) {
              Real* const Mrow = &M[npc2*(npr*ri + rj)];
              const Real fac2 = fac1*rvy[ri]*rvx[rj];
              for (Int ci = 0; ci < npc; ++ci)
                for (Int cj = 0; cj < npc; ++cj)
                  Mrow[npc*ci + cj] += fac2*cvy[ci]*cvx[cj];
            }
        }
      }
    }
  }

  return true;
}

} // namespace slmm
