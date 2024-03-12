#include "slmm_vis.hpp"
#include "slmm_io.hpp"
#include "slmm_util.hpp"

#include "siqk_search.hpp"
#include "siqk_geometry.hpp"

#include <fstream>

namespace slmm {
namespace vis {

Reconstruction::Enum Reconstruction::convert (const std::string& s) {
  if (s == "constant") return Enum::constant;
  return Enum::bilin;
}

std::string Reconstruction::convert (const Enum e) {
  if (e == Enum::constant) return "constant";
  return "bilin";
}


struct SourceType {
  enum Enum { gll, physgrid };
};

typedef siqk::Octree<siqk::SphereGeometry, 10> Octree;

struct SearchFunctor {
  static constexpr Real pad = 1e-5;
  const AVec3s ps, nmls;
  const AIdxs c2n, c2nml;
  const Real* const p;
  Real bb[6];
  int ci_hit;

  SearchFunctor (const AVec3s& ps_, const AIdxs& c2n_,
                 const AVec3s& nml_, const AIdxs& c2nml_,
                 const Real p_[3])
    : ps(ps_), nmls(nml_), c2n(c2n_), c2nml(c2nml_), p(p_)
  {
    for (int j = 0; j < 3; ++j) bb[j  ] = p[j] - pad;
    for (int j = 0; j < 3; ++j) bb[j+3] = p[j] + pad;
    ci_hit = -1;
  }

  void operator() (const Int ci) {
    // Asssure there's always a hit in case numerics makes 'inside' return false
    // for all potential hits.
    if (ci_hit == -1) ci_hit = ci;
    const auto e = slice(c2n, ci);
    const auto en = slice(c2nml, ci);
    for (int i = 0; i < 4; ++i) {
      const Real* const a = slice(ps, e[i]);
      const Real* const nml = slice(nmls, en[i]);
      if ( ! siqk::SphereGeometry::inside(p, a, nml)) return;
    }
    ci_hit = ci;
  }
};

struct SparseMat {
  std::vector<Int> r, c;
  std::vector<Real> s;
};

// y = A*x
static void csr_mvp (const SparseMat& A, const Real* x, Real* y) {
  const Int nr = (Int) A.r.size() - 1;
# pragma omp parallel for
  for (Int r = 0; r < nr; ++r) {
    if (A.r[r+1] == A.r[r]) {
      y[r] = 0;
      continue;
    }
    Real accum = 0;
    Real min_x = 1e10, max_x = -1e10;
    for (Int j = A.r[r]; j < A.r[r+1]; ++j) {
      const auto xj = x[A.c[j]];
      accum += A.s[j] * xj;
      min_x = std::min(min_x, xj);
      max_x = std::max(max_x, xj);
    }
    // Bilinear interp is monotone but it fail to be so numerically. It's
    // important not to introduce new extrema in the output in this
    // diagnostic. So limit it.
    accum = std::min(max_x, std::max(accum, min_x));
    y[r] = accum;
  }
}

static void cclinspace (const Int n, const Real lo, const Real hi,
                        std::vector<Real>& x) {
  x.resize(n);
  for (Int i = 0; i < n; ++i) {
    const Real a = (i + 0.5)/n;
    x[i] = (1 - a)*lo + a*hi;
  }
}

static void
fill_row (const AVec3s& p, const AIdxs& c2n, const Int ci, const Int r,
          const Real* s, const SourceType::Enum src_type,
          const Reconstruction::Enum recon, SparseMat& op) {
  const auto& cell = slice(c2n, ci);
  Real* ls = &op.s[op.r[r]];
  if (src_type == SourceType::gll) {
    // Compute reference coordinates.
    Real a, b;
    siqk::sqr::Info info;
    siqk::sqr::calc_sphere_to_ref(p, cell, s, a, b, &info);
    // Fill matrix row.
    for (Int i = 0; i < 4; ++i) op.c[op.r[r]+i] = cell[i];
    if (recon == Reconstruction::Enum::bilin) {
      static const Real f = 0.25;
      ls[0] = f*(1-a)*(1-b); ls[1] = f*(1+a)*(1-b);
      ls[2] = f*(1+a)*(1+b); ls[3] = f*(1-a)*(1+b);
    } else {
      ls[0] = ls[1] = ls[2] = ls[3] = 0.25;
    }
  } else {
    op.c[op.r[r]] = ci;
    *ls = 1;
  }
}

// Map for const or bilin interp for lat-lon tensor-grid output. c2n should be
// the I/O cell->node map, i.e., subcells of the full element.
static void init_latlon_ioc2n_map (
  const AVec3s& p, const AIdxs& c2n, const Octree& ot,
  const std::vector<Real>& lats, const std::vector<Real>& lons,
  const SourceType::Enum src_type, const Reconstruction::Enum recon,
  SparseMat& op)
{
  // Compute edge normals for octree search.
  AVec3s nml;
  AIdxs c2nml;
  fill_normals(p, c2n, nml, c2nml);
  // Make a sparse CSR matrix implementing bilinear interpolation. We know most
  // of the structure of this matrix up front.
  const Int nlon = (Int) lons.size(), nr = (Int) lats.size() * nlon;
  const Int nnz_per_row = src_type == SourceType::gll ? 4 : 1;
  op.r.resize(nr+1);
  op.c.resize(nnz_per_row*nr);
  op.s.resize(nnz_per_row*nr);
  op.r[0] = 0;
  for (Int r = 0; r < nr; ++r) op.r[r+1] = op.r[r] + nnz_per_row;
# pragma omp parallel for
  for (Int r = 0; r < nr; ++r) {
    // Find the cell containing this point.
    Real s[3];
    slmm::ll2xyz(lats[r / nlon], lons[r % nlon], s[0], s[1], s[2]);
    SearchFunctor sf(p, c2n, nml, c2nml, s);
    ot.apply(sf.bb, sf);
    const Int ci = sf.ci_hit;
    assert(ci >= 0 && ci < nslices(c2n));
    fill_row(p, c2n, ci, r, s, src_type, recon, op);
  }
}

// Map for const or bilin interp for orthographic output. c2n should be the I/O
// cell->node map, i.e., subcells of the full element.
static void init_orthographic_ioc2n_map (
  const AVec3s& p, const AIdxs& c2n, const Octree& ot,
  const std::vector<Real>& xs, const std::vector<Real>& ys,
  const Real xhat[3], const Real yhat[3], const Real zhat[3],
  const Reconstruction::Enum recon, SparseMat& op)
{
  using geo = siqk::SphereGeometry;
  // Compute edge normals for octree search.
  AVec3s nml;
  AIdxs c2nml;
  fill_normals(p, c2n, nml, c2nml);
  // Make a sparse CSR matrix implementing bilinear interpolation. We know most
  // of the structure of this matrix up front.
  const Int nx = (Int) xs.size(), nr = (Int) ys.size() * nx;
  const Int nnz_per_row = 4;
  op.r.resize(nr+1);
  op.r[0] = 0;
  for (Int r = 0; r < nr; ++r) {
    const Int ix = r % nx, iy = r / nx;
    const Real xy_mag2 = square(xs[ix]) + square(ys[iy]);
    op.r[r+1] = op.r[r] + (xy_mag2 < 1 ? nnz_per_row : 0);
  }
  op.c.resize(op.r[nr], 0);
  op.s.resize(op.r[nr], 0);
# pragma omp parallel for
  for (Int r = 0; r < nr; ++r) {
    // Find the cell containing this point.
    const Int ix = r % nx, iy = r / nx;
    const Real xy_mag2 = square(xs[ix]) + square(ys[iy]);
    Real s[3] = {0};
    if (xy_mag2 < 1) {
      const Real z = std::sqrt(1 - xy_mag2);
      geo::axpy(xs[ix], xhat, s);
      geo::axpy(ys[iy], yhat, s);
      geo::axpy(z, zhat, s);
    } else {
      continue;
    }
    SearchFunctor sf(p, c2n, nml, c2nml, s);
    ot.apply(sf.bb, sf);
    const Int ci = sf.ci_hit;
    assert(ci >= 0 && ci < nslices(c2n));
    fill_row(p, c2n, ci, r, s, SourceType::gll, recon, op);
  }
}

struct MapToLatLon::Impl {
  std::vector<Real> lats, lons;
  Octree ot;
  SparseMat op;

  Impl (const AVec3s& cgll_p, const AIdxs& cgll_io_c2n,
        Int nlat, Int nlon, const SourceType::Enum src_type,
        const Reconstruction::Enum recon) {
    cclinspace(nlat, -M_PI/2, M_PI/2, lats);
    cclinspace(nlon, -M_PI, M_PI, lons);
    ot.init(a2ConstVec3s(cgll_p), a2ConstIdxs(cgll_io_c2n));
    init_latlon_ioc2n_map(cgll_p, cgll_io_c2n, ot, lats, lons, src_type, recon,
                          op);
  }

  // Input is CGLL data, e.g., from D2Cer::d2c. Output is lon-lat data in a
  // rectangle, with longitude the faster index.
  void remap (const Real* field_cgll, Real* field_ll) const {
    csr_mvp(op, field_cgll, field_ll);
  }
};

const std::vector<Real>& MapToLatLon::get_lons () const { return p->lons; }
const std::vector<Real>& MapToLatLon::get_lats () const { return p->lats; }
size_t MapToLatLon::get_x_size () const { return get_lons().size(); }
size_t MapToLatLon::get_y_size () const { return get_lats().size(); }

BilinGLLToLatLon
::BilinGLLToLatLon (const AVec3s& cgll_p, const AIdxs& cgll_io_c2n,
                    Int nlat, Int nlon, Reconstruction::Enum recon) {
  p = std::make_shared<Impl>(cgll_p, cgll_io_c2n, nlat, nlon, SourceType::gll,
                             recon);
}

PhysgridToLatLon
::PhysgridToLatLon (const AVec3s& cgll_p, const AIdxs& cgll_io_c2n,
                    Int nlat, Int nlon) {
  p = std::make_shared<Impl>(cgll_p, cgll_io_c2n, nlat, nlon, SourceType::physgrid,
                             Reconstruction::Enum::constant);
}

// Input is CGLL data, e.g., from D2Cer::d2c. Output is lon-lat data in a
// rectangle, with longitude the faster index.
void BilinGLLToLatLon::remap (const Real* field_cgll, Real* field_ll) const {
  p->remap(field_cgll, field_ll);
}

void PhysgridToLatLon::remap (const Real* field_pg, Real* field_ll) const {
  p->remap(field_pg, field_ll);
}

struct BilinGLLToOrthographic::Impl {
  Real zhat[3];
  std::vector<Real> xs, ys;
  Octree ot;
  SparseMat op;

  Impl (const AVec3s& cgll_p, const AIdxs& cgll_io_c2n,
        const Real xhat[3], const Real yhat[3],
        Int nx, Int ny, Reconstruction::Enum recon)
  {
    using geo = siqk::SphereGeometry;
    assert(std::abs(geo::norm2(xhat) - 1) < 1e-10);
    assert(std::abs(geo::norm2(yhat) - 1) < 1e-10);
    geo::cross(xhat, yhat, zhat);
    assert(std::abs(geo::norm2(zhat) - 1) < 1e-10);
    const auto r = 1;
    cclinspace(nx, -r, r, xs);
    cclinspace(ny, -r, r, ys);
    ot.init(a2ConstVec3s(cgll_p), a2ConstIdxs(cgll_io_c2n));
    init_orthographic_ioc2n_map(cgll_p, cgll_io_c2n, ot, xs, ys,
                                xhat, yhat, zhat, recon, op);
  }

  void remap (const Real* field_cgll, Real* field_xy) const {
    csr_mvp(op, field_cgll, field_xy);
  }
};

BilinGLLToOrthographic
::BilinGLLToOrthographic (const AVec3s& cgll_p, const AIdxs& cgll_io_c2n,
                          const Real xhat[3], const Real yhat[3],
                          Int nx, Int ny, Reconstruction::Enum recon)
{
  p = std::make_shared<Impl>(cgll_p, cgll_io_c2n, xhat, yhat, nx, ny, recon);
}

void BilinGLLToOrthographic::remap (const Real* field_cgll, Real* field_xy) const {
  p->remap(field_cgll, field_xy);
}

size_t BilinGLLToOrthographic::get_x_size () const { return p->xs.size(); }
size_t BilinGLLToOrthographic::get_y_size () const { return p->ys.size(); }

struct VisWriter::Impl {
  MapToArray::Ptr op;
  io::InternalWriter::Ptr iw;
  std::vector<Real> wrk;
  Impl (const MapToArray::Ptr& op_, const std::string& filename)
    : op(op_), iw(std::make_shared<io::InternalWriter>(filename)),
      wrk(op->get_x_size() * op->get_y_size())
  {}
};

VisWriter::VisWriter (const MapToArray::Ptr& op, const std::string& filename) {
  p = std::make_shared<Impl>(op, filename);
}

void VisWriter::write (const Real* field_cgll) {
  p->op->remap(field_cgll, p->wrk.data());
  const Int dims[] = {(Int) p->op->get_y_size(),
                      (Int) p->op->get_x_size()};
  p->iw->write_array(2, dims, p->wrk.data());
}

MapToArray::Ptr VisWriter::get_map () { return p->op; }

} // namespace vis
} // namespace slmm
