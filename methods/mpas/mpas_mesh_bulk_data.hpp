#ifndef INCLUDE_MPAS_MESH_BULK_DATA_HPP
#define INCLUDE_MPAS_MESH_BULK_DATA_HPP
#ifdef MPAS_HAVE_NETCDF

#include "mpas_mesh.hpp"
#include "mpas_distribution.hpp"

struct MpasMeshBulkData {
  typedef std::shared_ptr<MpasMeshBulkData> Ptr;

  MpasMeshBulkData();
  bool is_ok() const;
  MpasMesh::Ptr get_mpas_mesh();

  static MpasMeshBulkData::Ptr create_and_update_distribution(
    const std::string& fname, Distribution& d, const Int nlevel = -1);

private:
  struct Impl;
  std::shared_ptr<Impl> p;

  // Create MpasMeshBulkData structure, and update the Distribution structure
  // with additional edge data.
  friend MpasMeshBulkData::Ptr create_and_update_distribution(
    const std::string& fname, Distribution& d, const Int nlevel);
};

#else // MPAS_HAVE_NETCDF

struct MpasMeshBulkData {
  typedef std::shared_ptr<MpasMeshBulkData> Ptr;
  bool is_ok() const { return false; }
  MpasMesh::Ptr get_mpas_mesh () { return nullptr; }
  static MpasMeshBulkData::Ptr
  create_and_update_distribution (const std::string&, Distribution&, const Int {
    printf("MpasMeshBulkData not supported. Provide NETCDF when building.\n");
    return nullptr;
  }
};

#endif // MPAS_HAVE_NETCDF
#endif
