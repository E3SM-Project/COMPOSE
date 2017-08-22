#ifndef INCLUDE_QLT_HPP
#define INCLUDE_QLT_HPP

#include <mpi.h>

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include <map>

#include <Kokkos_Core.hpp>
#include "qlt_kokkos.hpp"

// QLT: Quasi-local tree-based non-iterative tracer density reconstructor for
//      mass conservation, shape preservation, and tracer consistency.
namespace qlt {
typedef int Int;
typedef size_t Size;
typedef double Real;

namespace impl { class NodeSets; }

class Parallel {
  MPI_Comm comm_;
public:
  typedef std::shared_ptr<Parallel> Ptr;
  Parallel(MPI_Comm comm) : comm_(comm) {}
  MPI_Comm comm () const { return comm_; }
  Int size() const;
  Int rank() const;
  Int root () const { return 0; }
  bool amroot () const { return rank() == root(); }
};

Parallel::Ptr make_parallel(MPI_Comm comm);

namespace tree {
// The caller builds a tree of these nodes to pass to QLT.
struct Node {
  typedef std::shared_ptr<Node> Ptr;
  const Node* parent; // (Can't be a shared_ptr: would be a circular dependency.)
  Int rank;           // Owning rank.
  Int cellidx;        // If a leaf, the cell to which this node corresponds.
  Int nkids;          // 0 at leaf, 1 or 2 otherwise.
  Node::Ptr kids[2];
  void* reserved;     // For internal use.
  Node () : parent(nullptr), rank(-1), cellidx(-1), nkids(0), reserved(nullptr) {}
};
} // namespace tree

template <typename ExeSpace>
class QLT {
public:
  typedef typename impl::DeviceType<ExeSpace>::type Device;
  typedef QLT<ExeSpace> Me;
  typedef std::shared_ptr<Me> Ptr;
  
  struct ProblemType {
    enum : Int { conserve = 1, shapepreserve = 1 << 1, consistent = 1 << 2 };
  };

  // Set up QLT topology and communication data structures based on a tree.
  QLT(const Parallel::Ptr& p, const Int& ncells, const tree::Node::Ptr& tree);

  void print(std::ostream& os) const;

  // Number of cells owned by this rank.
  Int nlclcells() const;

  // Cells owned by this rank, in order of local numbering. Thus,
  // gci2lci(gcis[i]) == i. Ideally, the caller never actually calls gci2lci(),
  // and instead uses the information from get_owned_glblcells to determine
  // local cell indices.
  void get_owned_glblcells(std::vector<Int>& gcis) const;

  // For global cell index cellidx, i.e., the globally unique ordinal associated
  // with a cell in the caller's tree, return this rank's local index for
  // it. This is not an efficient operation.
  Int gci2lci(const Int& gci) const;

  // Set up QLT tracer metadata. Once end_tracer_declarations is called, it is
  // an error to call declare_tracer again. Call declare_tracer in order of the
  // tracer index in the caller's numbering.
  void declare_tracer(int problem_type);

  void end_tracer_declarations();

  int get_problem_type(const Int& tracer_idx) const;

  Int get_num_tracers() const;

  // set_{rho,Q}: Set cell values prior to running the QLT algorithm.
  //   set_rho must be called before set_Q.
  //   lclcellidx is gci2lci(cellidx).
  //   Notation:
  //     rho: Total density.
  //       Q: Tracer density.
  //       q: Tracer mixing ratio = Q/rho.
  //      *m: Mass corresponding to the density; results from an integral over a
  //          region, such as a cell.
  KOKKOS_INLINE_FUNCTION
  void set_rho(const Int& lclcellidx,
               // Current total mass in this cell.
               const Real& rhom);

  KOKKOS_INLINE_FUNCTION
  void set_Q(const Int& lclcellidx, const Int& tracer_idx,
             // Current tracer mass in this cell.
             const Real& Qm,
             // Minimum and maximum permitted tracer mass in this cell.
             const Real& Qm_min, const Real& Qm_max,
             // If mass conservation is requested, provide the previous Qm,
             // which will be summed to give the desired global mass.
             const Real Qm_prev = -1);

  // Run the QLT algorithm with the values set by set_{rho,Q}.
  void run();

  // Get a cell's tracer mass Qm after the QLT algorithm has run.
  KOKKOS_INLINE_FUNCTION
  Real get_Q(const Int& lclcellidx, const Int& tracer_idx);

private:
  typedef Kokkos::View<Int*, Kokkos::LayoutLeft, Device> IntList;
  typedef impl::Const<IntList> ConstIntList;
  typedef impl::ConstUnmanaged<IntList> ConstUnmanagedIntList;

  static void init(const std::string& name, IntList& d,
                   typename IntList::HostMirror& h, size_t n);

  struct MetaDataBuilder {
    typedef std::shared_ptr<MetaDataBuilder> Ptr;
    std::vector<int> trcr2prob;
  };

  struct MetaData {
    enum : Int { nprobtypes = 4 };

    template <typename IntListT>
    struct Arrays {
      // trcr2prob(i) is the ProblemType of tracer i.
      IntListT trcr2prob;
      // bidx2trcr(prob2trcrptr(i) : prob2trcrptr(i+1)-1) is the list of
      // tracers having ProblemType index i. bidx2trcr is the permutation
      // from the user's tracer index to the bulk data's ordering (bidx).
      Int prob2trcrptr[nprobtypes+1];
      IntListT bidx2trcr;
      // Inverse of bidx2trcr.
      IntListT trcr2bidx;
      // Points to the start of l2r bulk data for each problem type, within a
      // slot.
      Int prob2bl2r[nprobtypes + 1];
      // Point to the start of l2r bulk data for each tracer, within a slot.
      IntListT trcr2bl2r;
      // Same for r2l bulk data.
      Int prob2br2l[nprobtypes + 1];
      IntListT trcr2br2l;
    };

    static int get_problem_type(const int& idx);
    
    // icpc doesn't let us use problem_type_ here, even though it's constexpr.
    static int get_problem_type_idx(const int& mask);

    static int get_problem_type_l2r_bulk_size(const int& mask);

    static int get_problem_type_r2l_bulk_size(const int& mask);

    struct CPT {
      // We could make the l2r buffer smaller by one entry, Qm. However, the
      // l2r comm is more efficient if it's done with one buffer. Similarly,
      // we separate the r2l data into a separate buffer for packing and MPI
      // efficiency.
      //   There are 7 possible problems.
      //   The only problem not supported is conservation alone. It makes very
      // little sense to use QLT for conservation alone.
      //   The remaining 6 fall into 4 categories of details. These 4 categories
      // are traceked by QLT; which of the original 6 problems being solved is
      // not important.
      enum {
        // l2r: rhom, (Qm_min, Qm, Qm_max)*; l2r, r2l: Qm*
        s  = ProblemType::shapepreserve,
        st = ProblemType::shapepreserve | ProblemType::consistent,
        // l2r: rhom, (Qm_min, Qm, Qm_max, Qm_prev)*; l2r, r2l: Qm*
        cs  = ProblemType::conserve | s,
        cst = ProblemType::conserve | st,
        // l2r: rhom, (q_min, Qm, q_max)*; l2r, r2l: Qm*
        t = ProblemType::consistent,
        // l2r: rhom, (q_min, Qm, q_max, Qm_prev)*; l2r, r2l: Qm*
        ct = ProblemType::conserve | t
      };
    };

    Arrays<typename ConstUnmanagedIntList::HostMirror> a_h;
    Arrays<ConstUnmanagedIntList> a_d;

    void init(const MetaDataBuilder& mdb);

  private:
    static constexpr Int problem_type_[] = { CPT::st, CPT::cst, CPT::t, CPT::ct };
    Arrays<typename IntList::HostMirror> a_h_;
    Arrays<IntList> a_d_;
  };

  struct BulkData {
    typedef Kokkos::View<Real*, Kokkos::LayoutLeft, Device> RealList;
    typedef impl::Unmanaged<RealList> UnmanagedRealList;

    UnmanagedRealList l2r_data, r2l_data;

    void init(const MetaData& md, const Int& nslots);

  private:
    RealList l2r_data_, r2l_data_;
  };

private:
  void init(const Parallel::Ptr& p, const Int& ncells, const tree::Node::Ptr& tree);

  void init_ordinals();

  KOKKOS_INLINE_FUNCTION
  static void solve_node_problem(const Int problem_type,
                                 const Real& rhom, const Real* pd, const Real& Qm,
                                 const Real& rhom0, const Real* k0d, Real& Qm0,
                                 const Real& rhom1, const Real* k1d, Real& Qm1);

private:
  Parallel::Ptr p_;
  // Tree and communication topology.
  std::shared_ptr<const impl::NodeSets> ns_;
  // Globally unique cellidx -> rank-local index.
  std::map<Int,Int> gci2lci_;
  // Temporary to collect caller's tracer information prior to calling
  // end_tracer_declarations().
  typename MetaDataBuilder::Ptr mdb_;
  // Constructed in end_tracer_declarations().
  MetaData md_;
  BulkData bd_;
};

namespace test {
struct Input {
  bool unittest, write;
  Int ncells, ntracers, tracer_type, nrepeat;
  bool pseudorandom, verbose;
};

Int run(const Parallel::Ptr& p, const Input& in);
} // namespace test
} // namespace qlt

// These are the definitions that must be visible in the calling translation
// unit, unless Cuda relocatable device code is enabled.
#include "qlt_inline.hpp"

#endif
