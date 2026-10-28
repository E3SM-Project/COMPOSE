#ifndef INCLUDE_MPAS_NETCDF_HPP
#define INCLUDE_MPAS_NETCDF_HPP

#include "mpas_mesh.hpp"

#include <memory>
#include <string>

struct Comm;

struct NetcdfWriter {
  typedef std::shared_ptr<NetcdfWriter> Ptr;

  struct Options {
    bool double_precision;
    Options ()
      : double_precision(true)
    {}
  };

  struct Descriptor {
    std::string name;
  };
  typedef std::shared_ptr<const Descriptor> Desc;

  // Returns a descriptor used in write_field.
  Desc add_cc_lonlat_time_slice(const std::string& name, const int lvl,
                                const int qi = -1);
  Desc add_cc_vertlon_time_slice(const std::string& name, const Real latitude_radian,
                                 const int qi = -1);
  Desc add_cc_vert_greatcircle_time_slice(const std::string& name, const Vec axis,
                                          const int qi = -1);
  void end_definition();

  void advance_time_to(const double t);
  void write_field(const Desc& descriptor, const Acr2& rho);
  // If rho is null, write Q as is; else, write Q/rho.
  void write_field(const Desc& descriptor, const Acr2* rho, const Acr3& Q);

  static Ptr create(const std::shared_ptr<Comm>& c, const MpasMesh::ConstPtr& mm,
                    const std::string& fprefix, Options o = Options());

  static bool unittest(const std::shared_ptr<Comm>& c, const MpasMesh::ConstPtr& mm);

private:
  struct Impl;
  std::shared_ptr<Impl> p;
};

#endif
