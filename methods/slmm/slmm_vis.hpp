#ifndef INCLUDE_SLMM_VIS
#define INCLUDE_SLMM_VIS

#include <memory>

#include "slmm_defs.hpp"

namespace slmm {
namespace vis {

struct Reconstruction {
  enum class Enum { bilin, constant };
  static Enum convert(const std::string& s);
  static std::string convert(const Enum e);
};

struct MapToArray {
  typedef std::shared_ptr<MapToArray> Ptr;

  virtual void remap(const Real* field_from, Real* field_ll) const = 0;

  virtual size_t get_x_size() const = 0;
  virtual size_t get_y_size() const = 0;
};

class MapToLatLon : public MapToArray {
protected:
  class Impl;
  std::shared_ptr<Impl> p;
  
public:
  typedef std::shared_ptr<MapToLatLon> Ptr;

  const std::vector<Real>& get_lons() const;
  const std::vector<Real>& get_lats() const;

  size_t get_x_size() const override;
  size_t get_y_size() const override;
};

struct BilinGLLToLatLon : public MapToLatLon {
  typedef std::shared_ptr<BilinGLLToLatLon> Ptr;

  BilinGLLToLatLon(const AVec3s& cgll_p, const AIdxs& cgll_io_c2n,
                   Int nlat,  // cc(linspace(-pi/2, pi/2, nlat+1))
                   Int nlon,  // cc(linspace(-pi  , pi  , nlon+1))
                   Reconstruction::Enum recon = Reconstruction::Enum::bilin);

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

struct BilinGLLToOrthographic : public MapToArray {
  typedef std::shared_ptr<BilinGLLToOrthographic> Ptr;

  BilinGLLToOrthographic(const AVec3s& cgll_p, const AIdxs& cgll_io_c2n,
                         // Image frame.
                         const Real xhat[3], const Real yhat[3],
                         Int nx, Int ny,
                         Reconstruction::Enum recon = Reconstruction::Enum::bilin);

  // Input is CGLL data, e.g., from D2Cer::d2c. Output is xhat-yhat data in a
  // rectangle, with x the faster index.
  void remap(const Real* field_cgll, Real* field_xy) const override;

  size_t get_x_size() const override;
  size_t get_y_size() const override;

private:
  struct Impl;
  std::shared_ptr<Impl> p;
};

class VisWriter {
  class Impl;
  std::shared_ptr<Impl> p;

public:
  typedef std::shared_ptr<VisWriter> Ptr;

  VisWriter(const MapToArray::Ptr& bilin, const std::string& filename);
  void write(const Real* field_cgll);

  MapToArray::Ptr get_map();
};

} // namespace vis
} // namespace slmm

#endif
