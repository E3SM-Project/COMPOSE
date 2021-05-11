#ifndef INCLUDE_SLMM_IO_HPP
#define INCLUDE_SLMM_IO_HPP

#include "slmm_defs.hpp"

#include <vector>
#include <map>
#include <memory>

namespace netCDF {
class NcVar;
class NcFile;
}

namespace slmm {
namespace io {

struct FieldType { enum Enum { node, elem }; };

class NetcdfWriter {
  struct Field {
    std::string name;
    std::vector<std::shared_ptr<netCDF::NcVar> > ncvars;
    Field(const std::string& name, const Int dim);
  };

  Size nn_, nc_;
  Int time_idx_;
  double time_;
  bool define_done_;
  std::shared_ptr<netCDF::NcFile> ncf_;
  std::shared_ptr<netCDF::NcVar> time_v_;
  std::vector<Field> node_fields_, elem_fields_;
  typedef std::pair<FieldType::Enum, std::size_t> FieldIdx;
  std::map<std::string, FieldIdx> name2field_;

  void init(const AVec3s& p, const AIdxs& c2n,
            const std::string& out_fn, const Int np, const Int monotone_type);

public:
  typedef std::shared_ptr<NetcdfWriter> Ptr;

  // Open a Netcdf file for writing.
  NetcdfWriter(const AVec3s& p, const AIdxs& c2n,
               const std::string& out_fn,
               const Int np = 4, const Int monotone_type = 0);

  // Add and write a static field on the mesh to the file, such as nodal
  // weights.
  void add_and_write_nodal_static_field(const std::string& name,
                                        const double* field, const Int dim = 1);

  // Add time-dependent fields on the mesh to the file.
  void add_nodal_field(const std::string& name, const Int dim = 1);
  void add_element_field(const std::string& name, const Int dim = 1);

  // After adding all the time-dependent fields, end the definition phase,
  // providing the first time at which fields will be recorded.
  void end_definition();

  // Advance time forward before writing fields for a time step.
  void advance_time_to(const double t);

  // Write time-dependent field. If multidimensional, the fast index is the mesh
  // dimension.
  void write_field(const std::string& name, const double* field);
};

// vals must be preallocated.
void get_field_vals(const netCDF::NcFile& ncr, FieldType::Enum ft,
                    const int field_idx, const int time_idx, double* vals);

void get_field_names(
  const netCDF::NcFile& ncr, std::vector<std::string>& nodal_field_names,
  std::vector<std::string>& element_field_names);

Int get_np(const netCDF::NcFile& ncr);

class InternalWriter {
  class Impl;
  std::shared_ptr<Impl> p;

public:
  typedef std::shared_ptr<InternalWriter> Ptr;

  InternalWriter(const std::string& filename);
  void write_array(Int ndim, const Int* dims, const Real* data) const;
};

} // namespace io
} // namespace slmm

#endif
