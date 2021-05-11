#ifndef INCLUDE_SLMM_VIS
#define INCLUDE_SLMM_VIS

#include <memory>

#include "slmm_defs.hpp"

namespace slmm {
namespace vis {

class MapToLatLon {
protected:
  class Impl;
  std::shared_ptr<Impl> p;
  
public:
  typedef std::shared_ptr<MapToLatLon> Ptr;

  const std::vector<Real>& get_lons() const;
  const std::vector<Real>& get_lats() const;

  virtual void remap(const Real* field_from, Real* field_ll) const = 0;
};

struct BilinGLLToLatLon : public MapToLatLon {
  typedef std::shared_ptr<BilinGLLToLatLon> Ptr;

  BilinGLLToLatLon(const AVec3s& cgll_p, const AIdxs& cgll_io_c2n,
                   Int nlat,  // cc(linspace(-pi/2, pi/2, nlat+1))
                   Int nlon); // cc(linspace(-pi  , pi  , nlon+1))

  // Input is CGLL data, e.g., from D2Cer::d2c. Output is lon-lat data in a
  // rectangle, with longitude the faster index.
  void remap(const Real* field_cgll, Real* field_ll) const override;
};

struct PhysgridToLatLon : public MapToLatLon {
  typedef std::shared_ptr<PhysgridToLatLon> Ptr;

  PhysgridToLatLon(const AVec3s& cgll_p, const AIdxs& cgll_io_c2n,
                   Int nlat,  // cc(linspace(-pi/2, pi/2, nlat+1))
                   Int nlon); // cc(linspace(-pi  , pi  , nlon+1))

  // Input is PG data.
  void remap(const Real* field_pg, Real* field_ll) const override;
};

class VisWriter {
  class Impl;
  std::shared_ptr<Impl> p;

public:
  typedef std::shared_ptr<VisWriter> Ptr;

  VisWriter(const MapToLatLon::Ptr& bilin, const std::string& filename);
  void write(const Real* field_cgll);

  MapToLatLon::Ptr get_map();
};

} // namespace vis
} // namespace slmm

#endif
