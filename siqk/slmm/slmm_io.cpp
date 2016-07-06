#include "slmm_io.hpp"

#include <stdexcept>

#ifdef SLMM_HAVE_NETCDF
# include <netcdfcpp.h>
#endif

namespace slmm {
namespace io {

NetcdfWriter::NetcdfWriter (
  const Vec3s::HostMirror& p, const Idxs::HostMirror& c2n,
  const std::string& out_fn, const Int np, const Int monotone_type)
{
  init(p, c2n, out_fn, np, monotone_type);
}

void NetcdfWriter::init (
  const Vec3s::HostMirror& p, const Idxs::HostMirror& c2n,
  const std::string& out_fn, const Int np, const Int monotone_type)
{
#ifdef SLMM_HAVE_NETCDF
  nn_ = nslices(p);
  nc_ = nslices(c2n);

  time_idx_ = 0;
  time_ = 0;
  define_done_ = false;

  //todo Do I need this? NcError error(NcError::silent_nonfatal);
  ncf_ = std::make_shared<NcFile>(out_fn.c_str(), NcFile::Replace);
  if ( ! ncf_->is_valid())
    throw std::runtime_error(std::string("Could not open file ") + out_fn +
                             " for writing.");

  // Thank you, TempestRemap, for figuring out the Exodus stuff.
  static const int len_str = 33;
  auto nodes_dim = ncf_->add_dim("num_nodes", nn_);
  auto len_str_dim = ncf_->add_dim("len_string", len_str);
  auto time_dim = ncf_->add_dim("time_step");
  auto cells_dim = ncf_->add_dim("num_elem", nc_);
  auto num_el_blk_dim = ncf_->add_dim("num_el_blk", 1);
  auto nodes_per_cell_dim = ncf_->add_dim("num_nod_per_el1", szslice(c2n));
  auto att_block1_dim = ncf_->add_dim("num_att_in_blk1", 1);
  ncf_->add_dim("len_line", 81);
  ncf_->add_dim("num_dim", 3);
  ncf_->add_dim("num_el_in_blk1", nc_);
  ncf_->add_att("api_version", 4.98f);
  ncf_->add_att("version", 4.98f);
  ncf_->add_att("floating_point_word_size", 8);
  ncf_->add_att("file_size", 1);
  ncf_->add_att("title", "slmm::io::NetcdfWriter::init");

  ncf_->add_var("time_whole", ncDouble, time_dim);
  ncf_->add_var("eb_names", ncChar, num_el_blk_dim, len_str_dim);
  { // elem map
    std::vector<int> elem(nc_);
    for (Int i = 0; i < nc_; ++i) elem[i] = i+1;
    ncf_->add_var("elem_map", ncInt, cells_dim)->put(elem.data(), nc_);
  }
  { // c2n
    auto v = ncf_->add_var("connect1", ncInt, cells_dim, nodes_per_cell_dim);
    v->add_att("elem_type", "SHELL4");
    std::vector<int> connect(nc_*szslice(c2n));
    for (Int i = 0, k = 0; i < nslices(c2n); ++i)
      for (Int j = 0; j < szslice(c2n); ++j, ++k)
        connect[k] = c2n(i,j) + 1;
    v->set_cur(0, 0);
    v->put(connect.data(), nc_, szslice(c2n));
  }
  { // coords
    std::vector<double> buf(nn_);
    double* const d = buf.data();
    for (Int i = 0; i < nn_; ++i) d[i] = p(i,0);
    ncf_->add_var("coordx", ncDouble, nodes_dim)->put(d, nn_);
    for (Int i = 0; i < nn_; ++i) d[i] = p(i,1);
    ncf_->add_var("coordy", ncDouble, nodes_dim)->put(d, nn_);
    for (Int i = 0; i < nn_; ++i) d[i] = p(i,2);
    ncf_->add_var("coordz", ncDouble, nodes_dim)->put(d, nn_);
  }
  { // various other things
    int one = 1;
    ncf_->add_var("eb_status", ncInt, num_el_blk_dim)->put(&one, 1);
    auto v = ncf_->add_var("eb_prop1", ncInt, num_el_blk_dim);
    v->put(&one, 1);
    v->add_att("name", "ID");
    std::vector<double> buf(nc_, 1.0);
    v = ncf_->add_var("attrib1", ncDouble, cells_dim, att_block1_dim);
    v->put(buf.data(), nc_, 1);
  }

  add_att("np", np);
  add_att("monotone_type", monotone_type);
#else
  std::cerr << "Warning: NetcdfWriter::init: Netcdf was not compiled in.\n";
#endif
}

template <typename T>
void NetcdfWriter::add_att (const char* name, const T& val) {
#ifdef SLMM_HAVE_NETCDF
  ncf_->add_att(name, val);
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
  NcDim* const str_d = ncf_->get_dim("len_string");
  NcDim* const time_d = ncf_->get_dim("time_step");

  do {
    Int num_vars = 0;
    for (auto f: node_fields_)
      num_vars += static_cast<Int>(f.ncvars.size());
    if ( ! num_vars) break;

    NcDim* const nodes_d = ncf_->get_dim("num_nodes");
    NcDim* const nv_d = ncf_->add_dim("num_nod_var", num_vars);
    NcVar* const name_v = ncf_->add_var("name_nod_var", ncChar, nv_d, str_d);
    Int varno = 1;
    for (std::size_t i = 0; i < node_fields_.size(); ++i) {
      Field& f = node_fields_[i];
      if (f.ncvars.size() == 1) {
        name_v->set_cur(i, 0);
        name_v->put(f.name.c_str(), 1, f.name.size());

        std::stringstream ss;
        ss << "vals_nod_var" << varno++;
        f.ncvars[0] = ncf_->add_var(ss.str().c_str(), ncDouble, time_d, nodes_d);
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

    NcDim* const elem_d = ncf_->get_dim("num_elem");
    NcDim* const ev_d = ncf_->add_dim("num_elem_var", num_vars);
    NcVar* const name_v = ncf_->add_var("name_elem_var", ncChar, ev_d, str_d);
    Int varno = 1;
    for (std::size_t i = 0; i < elem_fields_.size(); ++i) {
      Field& f = elem_fields_[i];
      if (f.ncvars.size() == 1) {
        name_v->set_cur(i, 0);
        name_v->put(f.name.c_str(), 1, f.name.size());

        std::stringstream ss;
        ss << "vals_elem_var" << varno++ << "eb1";
        f.ncvars[0] = ncf_->add_var(ss.str().c_str(), ncDouble, time_d, elem_d);
      } else {
        //todo dim > 1
        throw std::runtime_error("dim > 1 not impl'ed.");
      }
    }    
  } while (0);

  time_ = -1;
  time_idx_ = -1;
  time_v_ = ncf_->get_var("time_whole");

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
    throw std::runtime_error("Invalid field.");
  Field& f = it->second.first == FieldType::node ?
    node_fields_[it->second.second] : elem_fields_[it->second.second];
  assert(f.ncvars.size() == 1); //todo dim > 1
  f.ncvars[0]->set_rec(time_idx_);
  f.ncvars[0]->put_rec(field);
  ncf_->sync();
#endif
}

void NetcdfWriter::advance_time_to (const double t) {
#ifdef SLMM_HAVE_NETCDF
  ++time_idx_;
  if (t <= time_)
    throw std::runtime_error("t must be > current time.");
  time_ = t;
  time_v_->set_rec(time_idx_);
  time_v_->put_rec(&time_);
#endif
}

NetcdfWriter::Field::Field (const std::string& name, const Int dim)
  : name(name), ncvars(dim, nullptr)
{}

void get_field_vals (const NcFile& ncr, FieldType::Enum ft, const int field_idx,
                     const int time_idx, double* vals) {
#ifdef SLMM_HAVE_NETCDF
  std::stringstream ss;
  int nvals;
  if (ft == FieldType::node) {
    ss << "vals_nod_var" << field_idx + 1;
    NcDim* const nodes_dim = ncr.get_dim("num_nodes");
    nvals = nodes_dim->size();
  } else {
    ss << "vals_elem_var" << field_idx + 1 << "eb1";
    NcDim* const cell_dim = ncr.get_dim("num_elem");
    nvals = cell_dim->size();
  }
  NcVar* const f_v = ncr.get_var(ss.str().c_str());
  f_v->set_cur(time_idx, 0);
  f_v->get(vals, 1, nvals);
#endif
}

void get_field_names (
  const NcFile& ncr, std::vector<std::string>& node_names,
  std::vector<std::string>& elem_names)
{
#ifdef SLMM_HAVE_NETCDF
  NcDim* const str_d = ncr.get_dim("len_string");
  std::vector<char> str(str_d->size());
  str.back() = '\0';
  do {
    NcDim* const nv_d = ncr.get_dim("num_nod_var");
    if ( ! nv_d) break;
    NcVar* const name_v = ncr.get_var("name_nod_var");
    for (int i = 0; i < nv_d->size(); ++i) {
      name_v->set_cur(i, 0);
      name_v->get(str.data(), 1, str.size());
      node_names.push_back(std::string(str.data()));
    }
  } while (0);
  do {
    NcDim* const ev_d = ncr.get_dim("num_elem_var");
    if ( ! ev_d) break;
    NcVar* const name_v = ncr.get_var("name_elem_var");
    for (int i = 0; i < ev_d->size(); ++i) {
      name_v->set_cur(i, 0);
      name_v->get(str.data(), 1, str.size());
      elem_names.push_back(std::string(str.data()));
    }
  } while (0);
#endif
}

#ifdef SLMM_HAVE_NETCDF
static NcValues* get_att_val (const NcFile& ncr, const char* name) {
  NcAtt* att;
  NcValues* vals;
  if ( ! (att = ncr.get_att(name)) ||
       ! (vals = att->values()))
    throw std::runtime_error(std::string("No attribute ") + name);
  delete att;
  return vals;
}
#endif

Int get_np (const NcFile& ncr) {
#ifdef SLMM_HAVE_NETCDF
  NcValues* vals = get_att_val(ncr, "np");
  const Int np = vals->as_int(0);
  delete vals;
  return np;
#else
  return 0;
#endif
}

} // namespace io
} // namespace slmm
