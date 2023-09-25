#ifndef INCLUDE_SLMMIR_SLMMIR_HPP
#define INCLUDE_SLMMIR_SLMMIR_HPP

#include "slmm_defs.hpp"
#include "slmm_util.hpp"
#include "slmm_spf.hpp"
using namespace slmm;

using PlaneGeo = siqk::PlaneGeometry;
using SphereGeo = siqk::SphereGeometry;

// Methods to implement discrete mass conservation.
struct Dmc {
  enum Enum {
    // Enforce conservation relative to quadrature of basis functions on the
    // sphere, using one equality constraint. Limiter QP feasibility appears to
    // be assured. Maintains OOA.
    eq_sphere,
    // Same, but w.r.t. Homme def of mass (integral of basis function on plane,
    // times Jacobian at node). QP feasibility seems to be assured because
    // source local mass is equal to target local mass. Seems to drop OOA by 1.
    eq_homme,
    // Enforce using facet transport, which in exact arithmetic is DMC. QP
    // feasibility seems to be assured.
    facet,
    // Use facet transport, but to get a couple of digits extra of DMC in FP,
    // use an equality constraint.
    eq_facet,
    // Enforce DMC using a single global equality constraint. The only case of
    // interest is the Homme mass, as I want to see if a single global
    // constraint will leave the OOA unaffected. I strongly suspect QP
    // feasibility is no longer assured because source and target local masses
    // are no longer necessarily equal.
    glbl_eq_homme,
    // Don't enforce DMC.
    none
  };
  static std::string convert (const Enum& e) {
    switch (e) {
    case Enum::eq_sphere: return "equality-sphere";
    case Enum::eq_homme: return "equality-homme";
    case Enum::facet: return "facet";
    case Enum::eq_facet: return "equality-facet";
    case Enum::glbl_eq_homme: return "global-equality-homme";
    case Enum::none:
    default: return "none";
    }
  }
  static Enum convert (const std::string& s) {
    if (eq(s, "es", "equality-sphere")) return Enum::eq_sphere;
    if (eq(s, "eh", "equality-homme")) return Enum::eq_homme;
    if (eq(s, "f", "facet")) return Enum::facet;
    if (eq(s, "ef", "equality-facet")) return Enum::eq_facet;
    if (eq(s, "geh", "global-equality-homme")) return Enum::glbl_eq_homme;
    return Enum::none;
  }
  static std::string get_inputs () {
    return "{equality-sphere, equality-homme, facet, equality-facet, "
      "global-equality-homme, none}";
  }
  static bool is_equality_constrained (const Enum& e) {
    return (e == eq_sphere || e == eq_homme || e == eq_facet ||
            e == glbl_eq_homme);
  }
  static bool is_locally_constrained (const Enum& e) {
    return (e == eq_sphere || e == eq_homme || e == eq_facet);
  }
  static bool is_globally_constrained (const Enum& e) {
    return (e == glbl_eq_homme);
  }
  static bool is_facet (const Enum& e) {
    return (e == facet || e == eq_facet);
  }
  static bool use_homme_mass (const Enum& e) {
    return (e == eq_homme || e == glbl_eq_homme || is_facet(e));
  }
};

// Various mesh types to study subcell meshes.
struct MeshType {
  enum Enum {
    // The usual cubed-sphere geometric mesh.
    geometric,
    // A mesh made by constructing a cubed-sphere geometric mesh, then refining
    // it into a GLL mesh.
    gllsubcell,
    // A mesh made by constructing a cubed-sphere geometric mesh, then refining
    // it by splitting the reference cells uniformly. Would it be better to
    // split the physical cells uniformly by arc length?
    runisubcell
  };
  static std::string convert (const Enum& e) {
    switch (e) {
    case Enum::geometric: return "geometric";
    case Enum::gllsubcell: return "gllsubcell";
    case Enum::runisubcell: return "runisubcell";
    default: throw std::runtime_error("MeshType: Not a valid enum.");
    }
  }
  static Enum convert (const std::string& s) {
    if (eq(s, "g", "geometric")) return Enum::geometric;
    if (eq(s, "gllsc", "gllsubcell")) return Enum::gllsubcell;
    if (eq(s, "runisc", "runisubcell")) return Enum::runisubcell;
    throw std::runtime_error(std::string("MeshType: Not a valid string: ") + s);
  }
  static bool is_subcell (const Enum& e) { return e != geometric; }
};

struct Method {
  enum Enum { ir, cdg, isl,
              pisl, /* stabilized interpolation SL (ISL) */
              pislu /* unstable ISL for testing and comparison */ };
  static std::string convert (const Enum& e) {
    switch (e) {
    case Enum::ir: return "ir";
    case Enum::cdg: return "cdg";
    case Enum::isl: return "isl";
    case Enum::pisl: return "pisl";
    case Enum::pislu: return "pislu";
    default: throw std::runtime_error("Method: Not a valid enum.");
    }
  }
  static Enum convert (const std::string& s) {
    if (eq(s, "ir")) return Enum::ir;
    if (eq(s, "cdg")) return Enum::cdg;
    if (eq(s, "isl")) return Enum::isl;
    if (eq(s, "pisl")) return Enum::pisl;
    if (eq(s, "pislu")) return Enum::pislu;
    // Aliases.
    if (eq(s, "csl")) return Enum::isl;
    if (eq(s, "pcsl")) return Enum::pisl;
    if (eq(s, "pcslu")) return Enum::pislu;
    throw std::runtime_error(std::string("Method: Not a valid string: ") + s);
  }
  static bool is_isl (const Enum& e) { return e == isl || e == pisl || e == pislu; }
  static bool is_pisl (const Enum& e) { return e == pisl || e == pislu; }
  static bool is_stabilized (const Enum& e) { return e != pislu; }
  static bool is_ir (const Enum& e) { return e == ir; }
  // Is cell-integrated.
  static bool is_ci (const Enum& e) { return ! is_isl(e); }
};

struct Filter {
  enum Enum { none = 0, qlt, qlt_pve, caas, mn2, caas_pve, caas_node };
  static std::string convert (const Enum& e) {
    switch (e) {
    case Enum::none: return "none";
    case Enum::qlt: return "qlt";
    case Enum::qlt_pve: return "qlt-pve";
    case Enum::caas: return "caas";
    case Enum::mn2: return "mn2";
    case Enum::caas_pve: return "caas-pve";
    case Enum::caas_node: return "caas-node";
    default: throw std::runtime_error("Filter: Not a valid enum.");
    }
  }
  static Enum convert (const std::string& s) {
    if (eq(s, "none")) return Enum::none;
    if (eq(s, "qlt")) return Enum::qlt;
    if (eq(s, "qlt-pve")) return Enum::qlt_pve;
    if (eq(s, "caas")) return Enum::caas;
    if (eq(s, "mn2")) return Enum::mn2;
    if (eq(s, "caas-pve")) return Enum::caas_pve;
    if (eq(s, "caas-node")) return Enum::caas_node;
    throw std::runtime_error(std::string("Filter: Not a valid string: ") + s);
  }
  static bool is_positive_only (const Enum& e) {
    return e == Enum::qlt_pve;
  }
  static bool is_nodal(const Enum& e) {
    return e == Enum::caas_node;
  }
  static spf::MassRedistributor::Method to_mrd (const Enum& e) {
    switch (e) {
    case Enum::qlt:
    case Enum::qlt_pve: return spf::MassRedistributor::qlt;
    case Enum::caas:
    case Enum::caas_pve: return spf::MassRedistributor::caas;
    case Enum::mn2: return spf::MassRedistributor::mn2;
    case Enum::caas_node:
      throw std::runtime_error("Filter::to_mrd: caas_node not supported.");
    default: throw std::runtime_error("Filter::to_mrd: Can't convert.");
    }
  }
};

struct Limiter {
  enum Enum { none, mn2, caas, caags, qlt };
  static std::string convert (const Enum& e) {
    switch (e) {
    case Enum::none: return "none";
    case Enum::mn2: return "mn2";
    case Enum::caas: return "caas";
    case Enum::caags: return "caags";
    case Enum::qlt: return "qlt";
    default: throw std::runtime_error("Limiter: Not a valid enum.");
    }
  }
  static Enum convert (const std::string& s) {
    if (eq(s, "none")) return Enum::none;
    if (eq(s, "mn2")) return Enum::mn2;
    if (eq(s, "caas")) return Enum::caas;
    if (eq(s, "caags")) return Enum::caags;
    if (eq(s, "qlt")) return Enum::qlt;
    throw std::runtime_error(std::string("Limiter: Not a valid string: ") + s);
  }
};

struct IOType {
  enum Enum { netcdf, internal };
  static std::string convert (const Enum& e) {
    switch (e) {
    case Enum::netcdf: return "netcdf";
    case Enum::internal: return "internal";
    default: throw std::runtime_error("IOType: Not a valid enum.");
    }
  }
  static Enum convert (const std::string& s) {
    if (eq(s, "netcdf")) return Enum::netcdf;
    if (eq(s, "internal")) return Enum::internal;
    throw std::runtime_error(std::string("IOType: Not a valid string: ") + s);
  }
};

#define puf(m) "(" << #m << " " << (m) << ")"
#define pu(m) << " " << puf(m)
template <typename T>
static void prarr (const std::string& name, const T* const v, const size_t n) {
  std::cerr << name << " = [";
  for (size_t i = 0; i < n; ++i) std::cerr << " " << v[i];
  std::cerr << "];\n";
}

template <typename T, typename Int>
inline void copy (T* const d, const T* const s, const Int& n) {
  memcpy(d, s, n*sizeof(T));
}

template <typename T>
inline void swap (T p[2]) {
  T tmp = p[0];
  p[0] = p[1];
  p[1] = tmp;
}

#if defined EXPENSIVE_CHECKS
# pragma message "EXPENSIVE CHECKS"
# define assert_expensive(m) assert(m)
#else
# define assert_expensive(m)
#endif

class Timer {
public:
  enum Op { ts_setup, ts, ts_integrate, ts_remap, ts_rest, ts_error,
            ts_remap_T, ts_remap_project, ts_remap_node_jac,
            ts_isl_interp,
            ts_remap_T_geometry, ts_remap_T_crs, ts_remap_T_fill,
            total, NTIMERS };
  static inline void init () {
#ifdef SLMM_TIME
    for (int i = 0; i < NTIMERS; ++i) et_[i] = 0;
#endif
  }
  static inline void start (const Op op) {
#ifdef SLMM_TIME
    gettimeofday(&t_start_[op], 0);
#endif
  }
  static inline void stop (const Op op) {
#ifdef SLMM_TIME
    timeval t2;
    gettimeofday(&t2, 0);
    const timeval& t1 = t_start_[op];
    static const double us = 1.0e6;
    et_[op] += (t2.tv_sec*us + t2.tv_usec - t1.tv_sec*us - t1.tv_usec)/us;
#endif
  }
# define tpr(op) do {                                                   \
    printf("%-20s %10.3e %10.1f\n", #op, et_[op], 100*et_[op]/tot);      \
  } while (0)
  static void print () {
#ifdef SLMM_TIME
    const double tot = et_[total];
    tpr(ts_setup); tpr(ts); tpr(ts_integrate); tpr(ts_remap);
    tpr(ts_remap_T); tpr(ts_remap_T_geometry); tpr(ts_remap_T_crs);
    tpr(ts_remap_T_fill); tpr(ts_remap_node_jac); tpr(ts_remap_project);
    tpr(ts_isl_interp);
    tpr(ts_rest); tpr(ts_error);
    printf("%-20s %10.3e %10.1f\n", "total", et_[total], 100.0);
#endif
  }
#undef tpr
private:
#ifdef SLMM_TIME
  static timeval t_start_[NTIMERS];
  static double et_[NTIMERS];
#endif
};

#endif
