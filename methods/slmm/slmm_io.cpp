#include "slmm_io.hpp"
#include "slmm_util.hpp"

#include <fstream>

#ifdef SLMM_HAVE_NETCDF
# include <netcdf>
#endif
using namespace netCDF;

#include <stdexcept>

namespace slmm {
namespace io {

NetcdfWriter::NetcdfWriter (
  const AVec3s& p, const AIdxs& c2n,
  const std::string& out_fn, const Int np, const Int monotone_type)
{
  init(p, c2n, out_fn, np, monotone_type);
}

void NetcdfWriter::init (
  const AVec3s& p, const AIdxs& c2n,
  const std::string& out_fn, const Int np, const Int monotone_type)
{
#ifdef SLMM_HAVE_NETCDF
  nn_ = nslices(p);
  nc_ = nslices(c2n);

  time_idx_ = 0;
  time_ = 0;
  define_done_ = false;

  try {
    ncf_ = std::make_shared<NcFile>(out_fn.c_str(), NcFile::replace);
  } catch (...) {
    throw std::runtime_error(std::string("Could not open file ") + out_fn +
                             " for writing.");
  }

  static const int len_str = 33;
  auto nodes_dim = ncf_->addDim("num_nodes", nn_);
  auto len_str_dim = ncf_->addDim("len_string", len_str);
  auto time_dim = ncf_->addDim("time_step");
  auto cells_dim = ncf_->addDim("num_elem", nc_);
  auto num_el_blk_dim = ncf_->addDim("num_el_blk", 1);
  auto nodes_per_cell_dim = ncf_->addDim("num_nod_per_el1", szslice(c2n));
  auto att_block1_dim = ncf_->addDim("num_att_in_blk1", 1);
  ncf_->addDim("len_line", 81);
  ncf_->addDim("num_dim", 3);
  ncf_->addDim("num_el_in_blk1", nc_);
  ncf_->putAtt("api_version", ncDouble, 4.98);
  ncf_->putAtt("version", ncDouble, 4.98);
  ncf_->putAtt("floating_point_word_size", ncInt, 8);
  ncf_->putAtt("file_size", ncInt, 1);
  ncf_->putAtt("title", "slmm::io::NetcdfWriter");

  ncf_->addVar("time_whole", ncDouble, time_dim);
  ncf_->addVar("eb_names", ncChar, {num_el_blk_dim, len_str_dim});
  { // elem map
    std::vector<int> elem(nc_);
    for (Int i = 0; i < nc_; ++i) elem[i] = i+1;
    ncf_->addVar("elem_map", ncInt, cells_dim).putVar(elem.data());
  }
  { // c2n
    auto v = ncf_->addVar("connect1", ncInt, {cells_dim, nodes_per_cell_dim});
    v.putAtt("elem_type", "SHELL4");
    std::vector<int> connect(nc_*szslice(c2n));
    for (Int i = 0, k = 0; i < nslices(c2n); ++i)
      for (Int j = 0; j < szslice(c2n); ++j, ++k)
        connect[k] = c2n(i,j) + 1;
    v.putVar(connect.data());
  }
  { // coords
    std::vector<double> buf(nn_);
    double* const d = buf.data();
    for (Int i = 0; i < nn_; ++i) d[i] = p(i,0);
    ncf_->addVar("coordx", ncDouble, nodes_dim).putVar(d);
    for (Int i = 0; i < nn_; ++i) d[i] = p(i,1);
    ncf_->addVar("coordy", ncDouble, nodes_dim).putVar(d);
    for (Int i = 0; i < nn_; ++i) d[i] = p(i,2);
    ncf_->addVar("coordz", ncDouble, nodes_dim).putVar(d);
  }
  { // various other things
    int one = 1;
    ncf_->addVar("eb_status", ncInt, num_el_blk_dim).putVar(&one);
    auto v = ncf_->addVar("eb_prop1", ncInt, num_el_blk_dim);
    v.putVar(&one);
    v.putAtt("name", "ID");
    std::vector<double> buf(nc_, 1.0);
    v = ncf_->addVar("attrib1", ncDouble, {cells_dim, att_block1_dim});
    v.putVar(buf.data());
  }

  ncf_->putAtt("np", ncInt, np);
  ncf_->putAtt("monotone_type", ncInt, monotone_type);
#else
  std::cerr << "Warning: NetcdfWriter::init: Netcdf was not compiled in.\n";
#endif
}

void NetcdfWriter::add_and_write_nodal_static_field (
  const std::string& name, const double* field, const Int dim)
{
#ifdef SLMM_HAVE_NETCDF
  if (dim != 1) throw std::runtime_error("dim != 1 is not impl'ed");
  auto nodes_dim = ncf_->getDim("num_nodes");
  ncf_->addVar(name.c_str(), ncDouble, nodes_dim).putVar(field);
#endif
}

void NetcdfWriter::add_nodal_field (const std::string& name, const Int dim) {
#ifdef SLMM_HAVE_NETCDF
  if (define_done_)
    throw std::runtime_error(
      "Can't add a new field after end_definition() was called.");
  const auto& it = name2field_.find(name);
  if (it != name2field_.end())
    throw std::runtime_error("Field name was already added.");
  name2field_[name] = FieldIdx(FieldType::node, node_fields_.size());
  node_fields_.push_back(Field(name, dim));
#endif
}

void NetcdfWriter::add_element_field (const std::string& name, const Int dim) {
#ifdef SLMM_HAVE_NETCDF
  if (define_done_)
    throw std::runtime_error(
      "Can't add a new field after end_definition() was called.");
  const auto& it = name2field_.find(name);
  if (it != name2field_.end())
    throw std::runtime_error("Field name was already added.");
  name2field_[name] = FieldIdx(FieldType::elem, elem_fields_.size());
  elem_fields_.push_back(Field(name, dim));
#endif
}

void NetcdfWriter::end_definition () {
#ifdef SLMM_HAVE_NETCDF
  NcDim str_d = ncf_->getDim("len_string");
  NcDim time_d = ncf_->getDim("time_step");

  do {
    Int num_vars = 0;
    for (auto f: node_fields_)
      num_vars += static_cast<Int>(f.ncvars.size());
    if ( ! num_vars) break;

    NcDim nodes_d = ncf_->getDim("num_nodes");
    NcDim nv_d = ncf_->addDim("num_nod_var", num_vars);
    NcVar name_v = ncf_->addVar("name_nod_var", ncChar, {nv_d, str_d});
    Int varno = 1;
    for (std::size_t i = 0; i < node_fields_.size(); ++i) {
      Field& f = node_fields_[i];
      if (f.ncvars.size() == 1) {
        name_v.putVar({i, 0}, {1, f.name.size()}, f.name.c_str());

        std::stringstream ss;
        ss << "vals_nod_var" << varno++;
        f.ncvars[0] = std::make_shared<NcVar>(
          ncf_->addVar(ss.str().c_str(), ncDouble, {time_d, nodes_d}));
      } else {
        //todo dim > 1
        throw std::runtime_error("dim > 1 not impl'ed.");
      }
    }
  } while (0);

  do {
    Int num_vars = 0;
    for (auto f: elem_fields_)
      num_vars += static_cast<Int>(f.ncvars.size());
    if ( ! num_vars) break;

    NcDim elem_d = ncf_->getDim("num_elem");
    NcDim ev_d = ncf_->addDim("num_elem_var", num_vars);
    NcVar name_v = ncf_->addVar("name_elem_var", ncChar, {ev_d, str_d});
    Int varno = 1;
    for (std::size_t i = 0; i < elem_fields_.size(); ++i) {
      Field& f = elem_fields_[i];
      if (f.ncvars.size() == 1) {
        name_v.putVar({i, 0}, {1, f.name.size()}, f.name.c_str());

        std::stringstream ss;
        ss << "vals_elem_var" << varno++ << "eb1";
        f.ncvars[0] = std::make_shared<NcVar>(
          ncf_->addVar(ss.str().c_str(), ncDouble, {time_d, elem_d}));
      } else {
        //todo dim > 1
        throw std::runtime_error("dim > 1 not impl'ed.");
      }
    }    
  } while (0);

  time_ = -1;
  time_idx_ = -1;
  time_v_ = std::make_shared<NcVar>(ncf_->getVar("time_whole"));

  define_done_ = true;
#endif
}

static void check_state (const Int time_idx, const bool define_done) {
#ifdef SLMM_HAVE_NETCDF
  if (time_idx == -1)
    throw std::runtime_error(
      "Need to advance_time_to before writing fields.");
  if ( ! define_done)
    throw std::runtime_error(
      "Can't write a field until end_definition() is called.");
#endif
}

void NetcdfWriter::write_field (const std::string& name, const double* field) {
#ifdef SLMM_HAVE_NETCDF
  check_state(time_idx_, define_done_);
  const auto& it = name2field_.find(name);
  if (it == name2field_.end())
    throw std::runtime_error("Invalid field " + name);
  const bool is_node = it->second.first == FieldType::node;
  Field& f = is_node ? node_fields_[it->second.second] :
    elem_fields_[it->second.second];
  assert(f.ncvars.size() == 1); //todo dim > 1
  f.ncvars[0]->putVar({(size_t) time_idx_, 0l},
                      {1l, (size_t) (is_node ? nn_ : nc_)},
                      field);
  ncf_->sync();
#endif
}

void NetcdfWriter::advance_time_to (const double t) {
#ifdef SLMM_HAVE_NETCDF
  ++time_idx_;
  if (t <= time_)
    throw std::runtime_error("t must be > current time.");
  time_ = t;
  time_v_->putVar({(size_t) time_idx_}, &time_);
#endif
}

NetcdfWriter::Field::Field (const std::string& name, const Int dim)
  : name(name), ncvars(dim)
{}

void get_field_vals (const NcFile& ncr, FieldType::Enum ft, const int field_idx,
                     const int time_idx, double* vals) {
#ifdef SLMM_HAVE_NETCDF
  std::stringstream ss;
  int nvals;
  if (ft == FieldType::node) {
    ss << "vals_nod_var" << field_idx + 1;
    NcDim nodes_dim = ncr.getDim("num_nodes");
    nvals = nodes_dim.getSize();
  } else {
    ss << "vals_elem_var" << field_idx + 1 << "eb1";
    NcDim cell_dim = ncr.getDim("num_elem");
    nvals = cell_dim.getSize();
  }
  NcVar f_v = ncr.getVar(ss.str().c_str());
  f_v.getVar({(size_t) time_idx, 0l}, {1l, (size_t) nvals}, vals);
#endif
}

void get_field_names (
  const NcFile& ncr, std::vector<std::string>& node_names,
  std::vector<std::string>& elem_names)
{
#ifdef SLMM_HAVE_NETCDF
  NcDim str_d = ncr.getDim("len_string");
  std::vector<char> str(str_d.getSize());
  str.back() = '\0';
  do {
    NcDim nv_d = ncr.getDim("num_nod_var");
    if (nv_d.isNull()) break;
    NcVar name_v = ncr.getVar("name_nod_var");
    for (size_t i = 0; i < nv_d.getSize(); ++i) {
      name_v.getVar({i, 0l}, {1l, str.size()}, str.data());
      node_names.push_back(std::string(str.data()));
    }
  } while (0);
  do {
    NcDim ev_d = ncr.getDim("num_elem_var");
    if (ev_d.isNull()) break;
    NcVar name_v = ncr.getVar("name_elem_var");
    for (size_t i = 0; i < ev_d.getSize(); ++i) {
      name_v.getVar({i, 0l}, {1l, str.size()}, str.data());
      elem_names.push_back(std::string(str.data()));
    }
  } while (0);
#endif
}

#ifdef SLMM_HAVE_NETCDF
template <typename T>
static T get_att_val (const NcFile& ncr, const char* name) {
  NcGroupAtt att = ncr.getAtt(name);
  if (att.isNull())
    throw std::runtime_error(std::string("No attribute ") + name);
  T value;
  att.getValues(&value);
  return value;
}
#endif

Int get_np (const NcFile& ncr) {
#ifdef SLMM_HAVE_NETCDF
  const Int np = get_att_val<Int>(ncr, "np");
  return np;
#else
  return 0;
#endif
}

struct InternalWriter::Impl {
  std::ofstream os;
  Impl (const std::string& filename)
    : os(filename.c_str(), std::ofstream::binary)
  {}
};

InternalWriter::InternalWriter (const std::string& filename) {
  p = std::make_shared<Impl>(filename);
}

void InternalWriter::write_array (Int ndim, const Int* dims, const Real* data) const {
  write(p->os, ndim);
  p->os.write((const char*) dims, ndim*sizeof(ndim));
  Int n = 1;
  for (Int i = 0; i < ndim; ++i) n *= dims[i];
  p->os.write((const char*) data, n*sizeof(Real));
}

} // namespace io
} // namespace slmm
