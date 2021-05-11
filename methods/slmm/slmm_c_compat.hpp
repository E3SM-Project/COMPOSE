#ifndef INCLUDE_SLMM_C_COMPAT_HPP
#define INCLUDE_SLMM_C_COMPAT_HPP

// Interface to some parts of SLMM for use in C or Fortran.

extern "C" {
  void init();
  void finalize();

  int solve_1eq_bc_qp_c(
    const int n, const double* w, const double* a, const double b,
    const double* xlo, const double* xhi, const double* y, double* x,
    const int max_its);

  void make_mesh(const int ne, const int np, const int nonuni, double* p, int* e);
  void make_latlon_quads(const int ncell, const double* p, const int* e, double* ll);
}

#endif
