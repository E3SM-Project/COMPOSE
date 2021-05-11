#ifndef INCLUDE_SLMMIR_SNAPSHOT_HPP
#define INCLUDE_SLMMIR_SNAPSHOT_HPP

#include "slmm_defs.hpp"
#include "slmm_util.hpp"
#include "slmm_gallery.hpp"
using namespace slmm;

#include "slmmir_mesh.hpp"

#include <memory>

struct Snapshot {
  typedef std::shared_ptr<Snapshot> Ptr;
  std::vector<gallery::InitialCondition::Shape> ics;
  std::vector<std::vector<Real> > q;
  Mesh::Ptr m;
};

void check(const Snapshot& s, const Mesh& m, const Real* q, const Real* w, const Int dnn);

#endif
