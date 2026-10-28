#ifndef INCLUDE_MPAS_DISTRIBUTION_HPP
#define INCLUDE_MPAS_DISTRIBUTION_HPP

#include "mpas_partition.hpp"
#include "mpas_mpi.hpp"

#include <map>

struct Distribution {
  typedef std::shared_ptr<Distribution> Ptr;
  typedef std::shared_ptr<const Distribution> ConstPtr;

  //// General data.
  mpi::Parallel::Ptr p;
  // Remotes with which we can comm.
  // pid(0) = my pid
  // pid(1:nremote-1) = pids of remotes
  Ai1 pids;
  Int szhalo; // number of cell layers in halo

  // Distribution data for an entity ('ent') type, cell or edge.
  struct Data {
    Int nown;  // number of own ents
    Int nhalo; // number of halo ents
    // pidptrs(0) = 0
    // pidptrs(1) = num owned ents
    // pidptrs(nremote) = num own + halo ents
    Ai1 pidptrs;
    // ent gids in my lid order; (l,g)ids are base-0 here
    Ai1 gids;
    // rlids(i) is the remote's lid for my lid nown+i
    Ai1 rlids; // (0:nchalo-1)

    // These additional data support full-halo exchanges for testing.
    // olids_ptrs(0:nremote-1) are pointers into olids. We also use these
    // pointers for the upwind MPI pattern, but we don't need olids in that
    // case.
    Ai1 olids_ptrs;
    // olids(olids_ptrs(i):olids_ptrs(i+1)-1) are my own-ent LIDs for ents to
    // send to remote pids(i).
    Ai1 olids;
  };

  Data c;
  Data e;
  // 0 if own, otherwise distance >= 1 from halo outer edge.
  Ac1 cdhalo;

  // c2r and cg are for the global grid. In the future, we could replace this
  // with a multi-round topography discovery routine. For now, that has low
  // priority.
  //   This function fill the cell data but not the edge data.
  static Distribution::Ptr
  create(const mpi::Parallel::Ptr& p, const std::vector<Int>& c2r,
         const CellGraph& cg, const Int szhalo);
};

// Initialize the edge data in d using c2e and preliminary edge LIDs. In
// addition, reset c2e, egids, and eg2l to be consistent. This routine does not
// use MPI. IDs are base-1.
void init_my_edge_data(const Int maxedge, Int* const c2e, const Int* c2ne,
                       Int* egids, std::map<Int,Int>& eg2l, Int& neown,
                       Distribution& d);
// Once init_my_edge_data has been called, exchange edge data among processes to
// set up comm data structures.
void exchange_edge_data(Distribution& d);

#endif
