#include <cassert>

#include <array>
#include <vector>
#include <limits>

#include "islet_tables.hpp"
#include "islet_util.hpp"
#include "islet_isl.hpp"
#include "islet_xnodes_metrics.hpp"
#include "islet_npx.hpp"

extern "C" {
  void dgemm_(const char* transa, const char* transb, const int* m,
              const int* n, const int* k, const double* alpha, const double* a,
              const int* lda, const double* b, const int* ldb,
              const double* beta, double* c, const int* ldc);
  void dpotrf_(const char* uplo, const int* n, double* a, const int* lda,
               int* info);
  void dpotrs_(const char* uplo, const int* n, const int* nrhs, const double* a,
               const int* lda, double* b, const int* ldb, int* info);
  void dtrsm_(const char* side, const char* uplo, const char* transa, const char* diag,
              const int* n, const int* nrhs, const double* alpha, const double* a,
              const int* lda, double* b, const int* ldb);
  void dtrtrs_(const char* uplo, const char* trans, const char* diag,
               const int* n, const int* nrhs, double* a, const int* lda,
               double* b, const int* ldb, int* info);
  void dgeqrf_(const int* m, const int* n, double* a, const int* lda,
               double* tau, double* wrk, int* iwrk, int* info);
  void dormqr_(const char* side, const char* trans,
               const int* m, const int* n, const int* k,
               double* a, const int* lda,
               double* tau, double* c, const int* ldc,
               double* wrk, const int* iwrk, int* info);
}

namespace islet {
// C = alpha op(A) op(B) + beta C
void dgemm (char transa, char transb, int m, int nrhs, int n, double alpha,
            const double* a, int lda, const double* b, int ldb, double beta,
            const double* c, int ldc) {
  dgemm_(&transa, &transb, &m, &nrhs, &n, &alpha, const_cast<double*>(a), &lda,
         const_cast<double*>(b), &ldb, &beta, const_cast<double*>(c), &ldc);
}

int dpotrf (char uplo, int n, double* a, int lda) {
  int info;
  dpotrf_(&uplo, &n, a, &lda, &info);
  return info;
}

int dpotrs (char uplo, int n, int nrhs, const double* a, int lda, double* bx,
            int ldb) {
  int info;
  dpotrs_(&uplo, &n, &nrhs, const_cast<double*>(a), &lda, bx, &ldb, &info);
  return info;
}

void dtrsm (char side, char uplo, char transa, char diag, int n, int nrhs,
            double alpha, const double* a, int lda, double* bx, int ldb) {
  dtrsm_(&side, &uplo, &transa, &diag, &n, &nrhs, &alpha,
         const_cast<double*>(a), &lda, bx, &ldb);
}

int dtrtrs (char uplo, char trans, char diag, int n, int nrhs,
            double* a, int lda, double* b, int ldb) {
  int info;
  dtrtrs_(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, &info);
  return info;
}

// tau[min(m,n)], wrk[>= n]
int dgeqrf (int m, int n, double* a, int lda,
            double* tau, double* wrk, int iwrk) {
  int info;
  dgeqrf_(&m, &n, a, &lda, tau, wrk, &iwrk, &info);
  return info;
}

// tau[min(m,n)], wrk[>= max(m,n)]
int dormqr (char side, char trans, int m, int n, int k, double* a, int lda,
            double* tau, double* c, int ldc, double* wrk, int iwrk) {
  int info;
  dormqr_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, wrk, &iwrk, &info);
  return info;
}

struct GllNatural : public Operator {
  virtual void eval (const Int& np, const Real& x, Real* const v) const override {
    eval_lagrange_poly(get_xnodes(np), np, x, v);
  }
};

struct GllOffsetNodalSubset : public Operator, public npxstab<Real> {
  virtual void eval (const Int& np, const Real& x, Real* const v) const override {
    npxstab<Real>::eval(np, x, v);
  }
};

void eval_offset_nodal_subset (
  const Int np, const Int nreg, const Int* subnp, const Int* os, const Real* xnodes,
  const Real& x, Real* const v)
{
  if (x > 0) {
    eval_offset_nodal_subset(np, nreg, subnp, os, xnodes, -x, v);
    for (int i = 0; i < np/2; ++i)
      std::swap(v[i], v[np-i-1]);
    return;
  }
  bool done = false;
  for (Int i = 0; i < nreg; ++i)
    if (x < xnodes[i+1]) {
      std::fill(v, v + np, 0);
      eval_lagrange_poly(xnodes + os[i], subnp[i], x, v + os[i]);
      done = true;
      break;
    }
  if ( ! done)
    eval_lagrange_poly(xnodes, np, x, v);
}

static void eval_offset (const Int& np, const Real* const xnodes,
                         const Int* const subnp, const Int* const offst,
                         const Real& x, Real* const v) {
  if (x > 0) {
    eval_offset(np, xnodes, subnp, offst, -x, v);
    for (int i = 0; i < np/2; ++i)
      std::swap(v[i], v[np-i-1]);
    return;
  }
  bool done = false;
  for (Int i = 0; i < np/2; ++i)
    if (x < xnodes[i+1]) {
      std::fill(v, v + np, 0);
      eval_lagrange_poly(xnodes + offst[i], subnp[i], x, v + offst[i]);
      done = true;
      break;
    }
  if ( ! done)
    eval_lagrange_poly(xnodes, np, x, v);
}

struct GllBest : public Operator {
  static void eval_np4 (const Real* const xnodes, const Real& x, Real* const y) {
    static const Real c1 = 0.306;
    if (x < xnodes[1] || x > xnodes[2]) {
      y[0] = y[3] = 0;
      const Int os = x < xnodes[1] ? 0 : 1;
      eval_lagrange_poly(xnodes + os, 3, x, y + os);
      Real y4[4];
      eval_lagrange_poly(xnodes, 4, x, y4);
      const Real x0 = 2*(1 - std::abs(x))/(1 - xnodes[2]) - 1;
      const Real a = (c1 + (0.5 - c1)*x0)*(x0 + 1);
      for (int i = 0; i < 4; ++i)
        y[i] = a*y[i] + (1 - a)*y4[i];
    } else
      eval_lagrange_poly(xnodes, 4, x, y);
  }

  virtual void eval (const Int& np, const Real& x, Real* const v) const override {
    const Real* xnodes = get_xnodes(np);
    switch (np) {
    case 4: eval_np4(xnodes, x, v); break; // 2
    case 5: { // 2
      const Int subnp[] = {3,4};
      const Int offst[] = {0,0};
      eval_offset(5, xnodes, subnp, offst, x, v);
    } break;
    case 6: { // 4
      const Int subnp[] = {5,5,6};
      const Int n0[] = { 0, 1, 2, 3, 4,  };
      const Int n1[] = { 0, 1, 2, 3,    5};
      const Int n2[] = { 0, 1, 2, 3, 4, 5};
      const Int* nodes[] = {n0,n1,n2};
      ::eval(6, true, xnodes, subnp, nodes, x, v);
    } break;
    case 7: { // 4
      const Int subnp[] = {5,5,6};
      const Int offst[] = {0,0,0};
      eval_offset(7, xnodes, subnp, offst, x, v);
    } break;
    case 8: { // 5
      const Int subnp[] = {6,6,7,6};
      const Int offst[] = {0,0,0,1};
      eval_offset(8, xnodes, subnp, offst, x, v);
    } break;
    case 9: { // 6
      const Int subnp[] = {7,8,8,7};
      const Int n0[] = { 0, 1, 2, 3, 4, 5,       8};
      const Int n1[] = { 0, 1, 2, 3, 4, 5,    7, 8};
      const Int n2[] = { 0, 1, 2, 3, 4, 5, 6,    8};
      const Int n3[] = {    1, 2, 3, 4, 5, 6, 7   };
      const Int* nodes[] = {n0,n1,n2,n3};
      ::eval(9, true, xnodes, subnp, nodes, x, v);
    } break;
    case 10: { // 6
      const Int subnp[] = {7,7,7,8,8};
      const Int offst[] = {0,0,0,0,1};
      eval_offset(10, xnodes, subnp, offst, x, v);
    } break;
    case 11: { // 7
      const Int subnp[] = {8,9,8,9,8};
      const Int offst[] = {0,0,0,0,1};
      eval_offset(11, xnodes, subnp, offst, x, v);
    } break;
    case 12: { // 8
      const Int subnp[] = {9,9,10,10,9,10};
      const Int offst[] = {0,0,0,0,1,1};
      eval_offset(12, xnodes, subnp, offst, x, v);
    } break;
    case 13: { // 9
      const Int subnp[] = {10,10,10,10,11,10};
      const Int offst[] = {0,0,0,0,0,1};
      eval_offset(13, xnodes, subnp, offst, x, v);
    } break;
    default: throw_if(true, "not impl'ed");
    }    
  }

  std::string get_basis_string (const Int& np) const override {
    switch (np) {
    case 5:  return "5 1 | 0 3: 0 1 2 | 1 4: 0 1 2 3";
    case 6:  return "6 1 | 0 5: 0 1 2 3 4 | 1 5: 0 1 2 3 5 | 2 6: 0 1 2 3 4 5";
    case 7:  return "7 1 | 0 5: 0 1 2 3 4 | 1 5: 0 1 2 3 4 | 2 6: 0 1 2 3 4 5";
    case 8:  return "8 1 | 0 6: 0 1 2 3 4 5 | 1 6: 0 1 2 3 4 5 | 2 7: 0 1 2 3 4 5 6 | 3 6: 1 2 3 4 5 6";
    case 9:  return "9 1 | 0 7: 0 1 2 3 4 5 8 | 1 8: 0 1 2 3 4 5 7 8 | 2 8: 0 1 2 3 4 5 6 8 | 3 7: 1 2 3 4 5 6 7";
    case 10: return "10 1 | 0 7: 0 1 2 3 4 5 6 | 1 7: 0 1 2 3 4 5 6 | 2 7: 0 1 2 3 4 5 6 | 3 8: 0 1 2 3 4 5 6 7 | 4 8: 1 2 3 4 5 6 7 8";
    case 11: return "11 1 | 0 8: 0 1 2 3 4 5 6 7 | 1 9: 0 1 2 3 4 5 6 7 8 | 2 8: 0 1 2 3 4 5 6 7 | 3 9: 0 1 2 3 4 5 6 7 8 | 4 8: 1 2 3 4 5 6 7 8";
    case 12: return "12 1 | 0 9: 0 1 2 3 4 5 6 7 8 | 1 9: 0 1 2 3 4 5 6 7 8 | 2 10: 0 1 2 3 4 5 6 7 8 9 | 3 10: 0 1 2 3 4 5 6 7 8 9 | 4 9: 1 2 3 4 5 6 7 8 9 | 5 10: 1 2 3 4 5 6 7 8 9 10";
    case 13: return "13 1 | 0 10: 0 1 2 3 4 5 6 7 8 9 | 1 10: 0 1 2 3 4 5 6 7 8 9 | 2 10: 0 1 2 3 4 5 6 7 8 9 | 3 10: 0 1 2 3 4 5 6 7 8 9 | 4 11: 0 1 2 3 4 5 6 7 8 9 10 | 5 10: 1 2 3 4 5 6 7 8 9 10";
    default: return "";
    }
  }
};

struct UniformOffsetNodalSubset : public Operator {
  virtual const Real* get_xnodes (const Int& np) const override {
    if (np < 2 || np > np_max+1) return nullptr;
    static Real xnode[np_max+1][np_max+1] = {0};
    if (xnode[np][0] == 0) {
      for (Int i = 0; i < np; ++i)
        xnode[np][i] = 2*(Real(i)/(np-1)) - 1;
    }
    return xnode[np];
  }

  virtual void eval (const Int& np, const Real& x, Real* const v) const override {
    const Real* xnodes = get_xnodes(np);
    switch (np) {
    case 2: {
      const Int subnp[] = {2};
      const Int offst[] = {0};
      eval_offset(2, xnodes, subnp, offst, x, v);
    } break;
    case 3: {
      const Int subnp[] = {3};
      const Int offst[] = {0};
      eval_offset(3, xnodes, subnp, offst, x, v);
    } break;
    case 4: {
      const Int subnp[] = {3,4};
      const Int offst[] = {0,0};
      eval_offset(4, xnodes, subnp, offst, x, v);
    } break;
    case 5: {
      const Int subnp[] = {3,4};
      const Int offst[] = {0,0};
      eval_offset(5, xnodes, subnp, offst, x, v);
    } break;
    case 6: {
      const Int subnp[] = {3,4,6};
      const Int offst[] = {0,0,0};
      eval_offset(6, xnodes, subnp, offst, x, v);
    } break;
    case 7: {
      const Int subnp[] = {3,4,4};
      const Int offst[] = {0,0,1};
      eval_offset(7, xnodes, subnp, offst, x, v);
    } break;
    case 8: {
      const Int subnp[] = {4,4,4,4};
      const Int offst[] = {0,0,1,2};
      eval_offset(8, xnodes, subnp, offst, x, v);
    } break;
    case 9: {
      const Int subnp[] = {4,4,4,4};
      const Int offst[] = {0,0,1,2};
      eval_offset(9, xnodes, subnp, offst, x, v);
    } break;
    case 10: {
      const Int subnp[] = {4,4,4,4,4};
      const Int offst[] = {0,0,1,2,3};
      eval_offset(10, xnodes, subnp, offst, x, v);
    } break;
    case 11: {
      const Int subnp[] = {4,4,4,4,4};
      const Int offst[] = {0,0,1,2,3};
      eval_offset(11, xnodes, subnp, offst, x, v);
    } break;
    case 12: {
      const Int subnp[] = {4,4,4,4,4,4};
      const Int offst[] = {0,0,1,2,3,4};
      eval_offset(12, xnodes, subnp, offst, x, v);
    } break;
    case 13: {
      const Int subnp[] = {4,4,4,4,4,4};
      const Int offst[] = {0,0,1,2,3,4};
      eval_offset(13, xnodes, subnp, offst, x, v);
    } break;
    default: throw_if(true, "not impl'ed");
    }    
  }
};

Operator::ConstPtr Operator::create (Operator::Method m) {
  switch (m) {
  case gll_natural: return std::make_shared<GllNatural>();
  case gll_offset_nodal_subset: return std::make_shared<GllOffsetNodalSubset>();
  case gll_best: return std::make_shared<GllBest>();
  case uniform_offset_nodal_subset: return std::make_shared<UniformOffsetNodalSubset>();
  default: throw_if(true, "Operator::create: not a method: " << m);
  }
  return nullptr;
}

Int unittest_eval () {
  Int nerr = 0;
  {
    GllOffsetNodalSubset o1;
    GllBest o2;
    for (const Int np : {5,7,8,10,11,12,13}) {
      const Int n = 100;
      Int ne = 0;
      for (Int i = 0; i <= n; ++i) {
        const Real x = 2*(Real(i)/n) - 1;
        Real v1[np_max], v2[np_max];
        o1.eval(np, x, v1);
        o2.eval(np, x, v2);
        for (Int j = 0; j < np; ++j) if (v1[j] != v2[j]) ++ne;
      }
      if (ne) printf("GllOffsetNodalSubset vs GllBest np %d failed\n", np);
      nerr += ne;
    }
  }
  return nerr;
}
} // namespace islet

using namespace islet;
extern "C" { // For python ctypes.
void get_xnodes (const Int method, const Int np, Real* xnodes) {
  const auto op = Operator::create(static_cast<Operator::Method>(method));
  const auto x = op->get_xnodes(np);
  for (Int i = 0; i < np; ++i) xnodes[i] = x[i];
}

void eval_interpolant (const Int method, const Int np, const Int nx,
                       // y is np x nx, np the fast index.
                       const Real* const x, Real* const y) {
  const auto op = Operator::create(static_cast<Operator::Method>(method));
  for (Int ix = 0; ix < nx; ++ix)
    op->eval(np, x[ix], y + np*ix);
}
} // extern "C"
