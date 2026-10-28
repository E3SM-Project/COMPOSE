#include "mpas_netcdf.hpp"
#include "mpas_comm.hpp"
#include "mpas_util.hpp"
#include "mpas_geometry.hpp"

#ifdef MPAS_HAVE_NETCDF
#include <vector>
#include <map>
#include <sstream>

#include <netcdf>
using namespace netCDF;

static void get_minmax_unit_z (const MpasMesh& mm, const int ci,
                               Real& z_min, Real& z_max) {
  z_min = 2; z_max = -2;
  const auto ne = mm.t.c2ne[ci];
  for (int i = 0; i < ne; ++i) {
    const auto vi = mm.t.get_c2v(ci,i);
    Real p[] = {mm.v.x[vi], mm.v.y[vi], mm.v.z[vi]};
    geo::normalize(p);
    z_min = std::min(z_min, p[2]);
    z_max = std::max(z_max, p[2]);
  }
}

// The plane goes through the origin and is given by the plane-normal 'nml'.
static void
get_minmax_extent_out_of_plane (const MpasMesh& mm, const Vec nml, const int ci,
                                Real& min_val, Real& max_val) {
  min_val = 2; max_val = -2;
  const auto ne = mm.t.c2ne[ci];
  for (int i = 0; i < ne; ++i) {
    const auto vi = mm.t.get_c2v(ci,i);
    Real p[] = {mm.v.x[vi], mm.v.y[vi], mm.v.z[vi]};
    geo::normalize(p);
    const Real val = geo::dot(nml, p);
    min_val = std::min(min_val, val);
    max_val = std::max(max_val, val);
  }
}

typedef NetcdfWriter::Desc Desc;

struct Field {
  enum Type : int { cc_lonlattime, cc_lonverttime, cc_gcverttime };
  // Superset of a field's description.
  Type type;
  int qi, lvl, lat_idx, axis_idx;
  NcVar var;
};

struct Latitude {
  Real latitude;
  NcDim dimension;
  std::vector<int> clids, cgids;
};

struct Axis {
  Vec axis;
  NcDim dimension;
  std::vector<int> clids, cgids;
};

struct NetcdfWriter::Impl {
  typedef std::shared_ptr<NcVar> NcVarPtr;

  std::shared_ptr<NcFile> ncf;
  // Time.
  Int time_idx;
  Real time;
  NcVar time_idx_v, time_v;
  // Descriptors.
  std::vector<Field> fields;
  std::map<Desc, size_t> desc2field;
  // Latitude data for vertical cross-sections given by latitude.
  std::vector<Latitude> lats;
  // Axis data for vertical cross-sections given by rotation axis.
  std::vector<Axis> axes;

  bool inited, def_done;
  Comm::Ptr c;
  MpasMesh::ConstPtr mm;

  void init (const Comm::Ptr& c_, const MpasMesh::ConstPtr& mm_,
             const std::string& fprefix, const NetcdfWriter::Options& o) {
    inited = false;
    def_done = false;
    c = c_;
    mm = mm_;

    try {
      std::stringstream ss;
      ss << fprefix << "_" << c->nrank() << "_" << c->rank() << ".nc";
      ncf = std::make_shared<NcFile>(ss.str().c_str(), NcFile::replace);
      inited = true;
    } catch (...) {
      return;
    }

    ncf->putAtt("num_mpi_rank", ncInt, c->nrank());
    ncf->putAtt("mpi_rank", ncInt, c->rank());
    ncf->putAtt("version", ncInt, 0);

    const auto time_dim = ncf->addDim("num_time");
    time_idx_v = ncf->addVar("time_idx", ncInt, time_dim);
    time_v = ncf->addVar("time", ncDouble, time_dim);
    time = -1;
    time_idx = -1;

    add_cc_metadata();
    add_vertical_metadata();
  }

  void add_cc_metadata() {
    const Int ncell = mm->ncown;
    const auto cell_dim = ncf->addDim("num_cell", ncell);
    auto cell_gid = ncf->addVar("cell_gid", ncInt, cell_dim);
    const Int* gids = get_cell_gids(*c);
    if (gids) {
      cell_gid.putVar(gids);
    } else {
      std::vector<Int> ids(ncell);
      for (Int i = 0; i < ncell; ++i) ids[i] = i;
      cell_gid.putVar(ids.data());
    }
  }

  void add_vertical_metadata () {
    const Int nmlvl = mm->l.nmax, nilvl = nmlvl+1;
    ncf->addDim("num_midpoint_level", nmlvl);
    ncf->addDim("num_interface_level", nilvl);
  }

  Desc create_Desc (const std::string& name, const Field::Type ftype, const int qi,
                    const int lvl, const int lat_idx, const int axis_idx) {
    auto d = std::make_shared<Descriptor>();
    std::stringstream s;
    s << name << "_ft" << int(ftype) << "_qi" << qi << "_lvl" << lvl
      << "_lati" << lat_idx << "_axisi" << axis_idx;
    d->name = s.str();
    return d;
  }

  int get_lat2idx (const Real lat) {
    for (size_t i = 0; i < lats.size(); ++i)
      if (lats[i].latitude == lat)
        return static_cast<int>(i);
    return -1;
  }

  int create_latitude (const Real lat) {
    lats.push_back(Latitude());
    const int idx = static_cast<int>(lats.size()) - 1;
    auto& l = lats.back();
    l.latitude = lat;
    fill_cc_lonverttime_arrays(lat, l.clids, l.cgids);
    std::stringstream s;
    s << "num_cgid_lati" << idx;
    l.dimension = ncf->addDim(s.str(), l.cgids.size());
    s.str("");
    s << "cell_gid_lati" << idx;
    auto cell_gid = ncf->addVar(s.str(), ncInt, l.dimension);
    cell_gid.putVar(l.cgids.data());
    return idx;
  }

  void fill_cc_lonverttime_arrays (const Real lat, std::vector<int>& clids,
                                   std::vector<int>& cgids) {
    clids.clear(); cgids.clear();
    const Real z0 = std::sin(lat);
    const Int* const own_cgids = get_cell_gids(*c);
    for (int ci = 0; ci < mm->ncown; ++ci) {
      Real min_z, max_z;
      get_minmax_unit_z(*mm, ci, min_z, max_z);
      if (z0 > max_z || z0 < min_z) continue;
      clids.push_back(ci);
      cgids.push_back(own_cgids[ci]);
    }
  }

  int get_axis2idx (const Vec axis) {
    for (size_t i = 0; i < axes.size(); ++i)
      if (geo::equal(axes[i].axis, axis))
        return static_cast<int>(i);
    return -1;
  }

  int create_axis (const Vec axis) {
    axes.push_back(Axis());
    const int idx = static_cast<int>(axes.size()) - 1;
    auto& a = axes.back();
    geo::copy(a.axis, axis);
    fill_cc_gcverttime_arrays(axis, a.clids, a.cgids);
    std::stringstream s;
    s << "num_cgid_axisi" << idx;
    a.dimension = ncf->addDim(s.str(), a.cgids.size());
    s.str("");
    s << "cell_gid_axisi" << idx;
    auto cell_gid = ncf->addVar(s.str(), ncInt, a.dimension);
    cell_gid.putVar(a.cgids.data());
    return idx;
  }

  void fill_cc_gcverttime_arrays (const Vec axis, std::vector<int>& clids,
                                  std::vector<int>& cgids) {
    clids.clear(); cgids.clear();
    const Int* const own_cgids = get_cell_gids(*c);
    for (int ci = 0; ci < mm->ncown; ++ci) {
      Real min_val, max_val;
      get_minmax_extent_out_of_plane(*mm, axis, ci, min_val, max_val);
      // Different signs mean the cell crosses the plane.
      if (min_val*max_val >= 0) continue;
      clids.push_back(ci);
      cgids.push_back(own_cgids[ci]);
    }
  }
};

NetcdfWriter::Desc NetcdfWriter
::add_cc_lonlat_time_slice (const std::string& name, const int lvl, const int qi) {
  const auto d = p->create_Desc(name, Field::cc_lonlattime, qi, lvl, -1, -1);
  p->desc2field[d] = p->fields.size();
  p->fields.push_back(Field());
  auto& f = p->fields.back();
  f.type = Field::Type::cc_lonlattime;
  f.qi = qi; f.lvl = lvl; f.lat_idx = -1;
  return d;
}

NetcdfWriter::Desc NetcdfWriter
::add_cc_vertlon_time_slice (const std::string& name, const Real latitude_radian,
                             const int qi) {
  int lat_idx = p->get_lat2idx(latitude_radian);
  if (lat_idx == -1) lat_idx = p->create_latitude(latitude_radian);
  const auto d = p->create_Desc(name, Field::cc_lonverttime, qi, -1, lat_idx, -1);
  p->desc2field[d] = p->fields.size();
  p->fields.push_back(Field());
  auto& f = p->fields.back();
  f.type = Field::Type::cc_lonverttime;
  f.qi = qi; f.lvl = -1; f.lat_idx = lat_idx; f.axis_idx = -1;
  return d;
}

NetcdfWriter::Desc NetcdfWriter
::add_cc_vert_greatcircle_time_slice (const std::string& name, const Vec axis_,
                                      const int qi) {
  Vec axis;
  geo::copy(axis, axis_);
  geo::normalize(axis);
  int axis_idx = p->get_axis2idx(axis);
  if (axis_idx == -1) axis_idx = p->create_axis(axis);
  const auto d = p->create_Desc(name, Field::cc_gcverttime, qi, -1, -1, axis_idx);
  p->desc2field[d] = p->fields.size();
  p->fields.push_back(Field());
  auto& f = p->fields.back();
  f.type = Field::Type::cc_gcverttime;
  f.qi = qi; f.lvl = -1; f.lat_idx = -1; f.axis_idx = axis_idx;
  return d;
}

void NetcdfWriter::end_definition () {
  { // Check that the descriptors are unique.
    bool ok = true;
    std::set<std::string> names;
    for (const auto& e : p->desc2field) {
      if (names.find(e.first->name) != names.end()) {
        printf("%s added multiple times.\n", e.first->name.c_str());
        ok = false;
      }
      names.insert(e.first->name);
    }
    throw_if( ! ok, "NetcdfWriter::end_definition: duplicate field names.");
  }
  auto& ncf = p->ncf;
  { // Add a variable for each array to be written.
    const auto cell_d = ncf->getDim("num_cell");
    const auto time_d = ncf->getDim("num_time");
    for (const auto& e : p->desc2field) {
      const auto& d = e.first;
      auto& f = p->fields[e.second];
      switch (f.type) {
      case Field::Type::cc_lonlattime:
        f.var = ncf->addVar(d->name.c_str(), ncDouble, {time_d, cell_d});
        break;
      case Field::Type::cc_lonverttime: {
        assert(f.lat_idx >= 0);
        const auto& lat = p->lats[f.lat_idx];
        const auto& vert_d = ncf->getDim("num_midpoint_level");
        f.var = ncf->addVar(d->name.c_str(), ncDouble,
                            {time_d, lat.dimension, vert_d});
      } break;
      case Field::Type::cc_gcverttime: {
        assert(f.axis_idx >= 0);
        const auto& a = p->axes[f.axis_idx];
        const auto& vert_d = ncf->getDim("num_midpoint_level");
        f.var = ncf->addVar(d->name.c_str(), ncDouble,
                            {time_d, a.dimension, vert_d});
      } break;
      }
    }
  }
  { // Add latitude array and metadata.
    const auto nlat = p->lats.size();
    std::vector<Real> lats(nlat);
    for (size_t i = 0; i < nlat; ++i)
      lats[i] = p->lats[i].latitude;
    const auto num_lat_d = ncf->addDim("num_lat", nlat);
    const auto lat_v = ncf->addVar("latitude", ncDouble, num_lat_d);
    lat_v.putVar(lats.data());
  }
  { // Add axis arrays and metadata.
    const auto nax = p->axes.size();
    std::vector<Real> axes(3*nax);
    for (size_t i = 0; i < nax; ++i)
      for (int d = 0; d < 3; ++d)
        axes[3*i+d] = p->axes[i].axis[d];
    const auto num_ax_d = ncf->addDim("num_axis", nax);
    const auto three_d = ncf->addDim("three", 3);
    const auto ax_v = ncf->addVar("axis", ncDouble, {num_ax_d, three_d});
    ax_v.putVar(axes.data());
  }
}

void NetcdfWriter::advance_time_to (const double t) {
  ++p->time_idx;
  throw_if(t <= p->time,
           "NetcdfWriter::advance_time_to: time must increase monotonically.");
  p->time = t;
  p->time_idx_v.putVar({size_t(p->time_idx)}, &p->time_idx);
  p->time_v.putVar({size_t(p->time_idx)}, &p->time);
}

static void
write_vert (const std::vector<int>& clids, const int nlvl, const int time_idx,
            const Acr2& rho, const Field& f) {
  const int ncell = static_cast<int>(clids.size());
  std::vector<Real> data(ncell*nlvl);
  for (int i = 0; i < ncell; ++i)
    for (int j = 0; j < nlvl; ++j)
      data[i*nlvl+j] = rho(j,clids[i]);
  f.var.putVar({size_t(time_idx), 0l, 0l},
               {1l, size_t(ncell), size_t(nlvl)},
               data.data());
}

void NetcdfWriter
::write_field (const Desc& descriptor, const Acr2& rho) {
  const auto it = p->desc2field.find(descriptor);
  throw_if(it == p->desc2field.end(),
           "NetcdfWriter::write_field: invalid descriptor "
           << descriptor->name);
  const auto& f = p->fields[it->second];
  switch (f.type) {
  case Field::Type::cc_lonlattime: {
    const auto ncown = p->mm->ncown;
    std::vector<Real> data(ncown);
    for (int i = 0; i < ncown; ++i)
      data[i] = rho(f.lvl,i);
    f.var.putVar({size_t(p->time_idx), 0l},
                 {1l, size_t(ncown)},
                 data.data());
  } break;
  case Field::Type::cc_lonverttime: {
    const auto& l = p->lats[f.lat_idx];
    write_vert(l.clids, p->mm->l.nmax, p->time_idx, rho, f);
  } break;
  case Field::Type::cc_gcverttime: {
    const auto& a = p->axes[f.axis_idx];
    write_vert(a.clids, p->mm->l.nmax, p->time_idx, rho, f);
  } break;
  }
  p->ncf->sync();
}

static void
write_vert (const std::vector<int>& clids, const int nlvl, const int time_idx,
            const Acr2* rhop, const Acr3& q, const Field& f) {
  const int ncell = static_cast<int>(clids.size());
  std::vector<Real> data(ncell*nlvl);
  if (rhop) {
    const auto& rho = *rhop;
    for (int i = 0; i < ncell; ++i)
      for (int j = 0; j < nlvl; ++j)
        data[i*nlvl+j] = q(f.qi,j,clids[i])/rho(j,clids[i]);
  } else {
    for (int i = 0; i < ncell; ++i)
      for (int j = 0; j < nlvl; ++j)
        data[i*nlvl+j] = q(f.qi,j,clids[i]);
  }
  f.var.putVar({size_t(time_idx), 0l, 0l},
               {1l, size_t(ncell), size_t(nlvl)},
               data.data());
}

void NetcdfWriter
::write_field (const Desc& descriptor, const Acr2* rhop, const Acr3& q) {
  const auto it = p->desc2field.find(descriptor);
  throw_if(it == p->desc2field.end(),
           "NetcdfWriter::write_field: invalid descriptor "
           << descriptor->name);
  const auto& f = p->fields[it->second];
  switch (f.type) {
  case Field::Type::cc_lonlattime: {
    const auto ncown = p->mm->ncown;
    std::vector<Real> data(ncown);
    if (rhop) {
      const auto& rho = *rhop;
      for (int i = 0; i < ncown; ++i)
        data[i] = q(f.qi,f.lvl,i)/rho(f.lvl,i);
    } else {
      for (int i = 0; i < ncown; ++i)
        data[i] = q(f.qi,f.lvl,i);
    }
    f.var.putVar({size_t(p->time_idx), 0l},
                 {1l, size_t(ncown)},
                 data.data());
  } break;
  case Field::Type::cc_lonverttime: {
    const auto& l = p->lats[f.lat_idx];
    write_vert(l.clids, p->mm->l.nmax, p->time_idx, rhop, q, f);
  } break;
  case Field::Type::cc_gcverttime: {
    const auto& a = p->axes[f.axis_idx];
    write_vert(a.clids, p->mm->l.nmax, p->time_idx, rhop, q, f);
  } break;
  }
  p->ncf->sync();
}

NetcdfWriter::Ptr NetcdfWriter
::create (const Comm::Ptr& c, const MpasMesh::ConstPtr& mm,
          const std::string& fprefix, Options o) {
  auto ncw = std::make_shared<NetcdfWriter>();
  ncw->p = std::make_shared<Impl>();
  ncw->p->init(c, mm, fprefix, o);
  return ncw->p->inited ? ncw : nullptr;
}

bool NetcdfWriter::unittest (const Comm::Ptr& c, const MpasMesh::ConstPtr& mm) {
  bool ok = true;
  auto ncw = NetcdfWriter::create(c, mm, "tmp/ncwtest");
  const int nq = 2, nlvl = 10, ncell = mm->t.ncell, ncown = mm->ncown;
  Ar3 q(nq, nlvl, ncell);
  const Int* cgids = get_cell_gids(*c);
  for (int i = 0; i < ncown; ++i)
    for (int k = 0; k < nlvl; ++k)
      for (int qi = 0; qi < nq; ++qi)
        q(qi,k,i) = cgids[i];
  const auto d0 = ncw->add_cc_lonlat_time_slice("q", 0, 0);
  const auto d1 = ncw->add_cc_vertlon_time_slice("q", 0, 0);
  const Real axis[] = {0,0,1};
  const auto d2 = ncw->add_cc_vert_greatcircle_time_slice("q", axis, 0);
  ncw->end_definition();
  ncw->advance_time_to(0);
  ncw->write_field(d0, nullptr, q);
  ncw->write_field(d1, nullptr, q);
  ncw->write_field(d2, nullptr, q);
  return ok;
}

#else
typedef NetcdfWriter::Desc Desc;
Desc NetcdfWriter
::add_cc_lonlat_time_slice (const std::string&, const int, const int)
{ return nullptr; }
Desc NetcdfWriter
::add_cc_vertlon_time_slice (const std::string&, const Real, const int)
{ return nullptr; }
Desc NetcdfWriter
::add_cc_vert_greatcircle_time_slice (const std::string&, const Vec, const int)
{ return nullptr; }
void NetcdfWriter::end_definition () {}
void NetcdfWriter::advance_time_to (const double) {}
void NetcdfWriter::write_field (const Desc&, const Acr2&) {}
NetcdfWriter::Ptr NetcdfWriter
::create (const Comm::Ptr&, const MpasMesh::ConstPtr&,
          const std::string&, Options) { return nullptr; }
bool NetcdfWriter
::unittest (const Comm::Ptr&, const MpasMesh::ConstPtr&) { return true; }
#endif
