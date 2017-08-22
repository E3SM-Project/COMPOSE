#ifndef INCLUDE_QLT_KOKKOS_HPP
#define INCLUDE_QLT_KOKKOS_HPP

namespace qlt {
namespace impl {
template <typename MemoryTraitsType, Kokkos::MemoryTraitsFlags flag>
using MemoryTraits = Kokkos::MemoryTraits<
  MemoryTraitsType::Unmanaged | MemoryTraitsType::RandomAccess |
  MemoryTraitsType::Atomic | flag>;

template <typename View>
using Unmanaged = Kokkos::View<
  typename View::data_type, typename View::array_layout,
  typename View::device_type, MemoryTraits<typename View::memory_traits,
                                           Kokkos::Unmanaged> >;
template <typename View>
using Const = Kokkos::View<
  typename View::const_data_type, typename View::array_layout,
  typename View::device_type, typename View::memory_traits>;
template <typename View>
using ConstUnmanaged = Const<Unmanaged<View> >;

template <typename ExeSpace>
struct DeviceType {
  typedef Kokkos::Device<typename ExeSpace::execution_space,
                         typename ExeSpace::memory_space> type;
};

#ifdef KOKKOS_HAVE_CUDA
typedef Kokkos::Device<Kokkos::CudaSpace::execution_space,
                       Kokkos::CudaSpace::memory_space> DefaultDeviceType;

template <> struct DeviceType<Kokkos::Cuda> {
  typedef DefaultDeviceType type;
};
#else
typedef Kokkos::Device<Kokkos::DefaultExecutionSpace::execution_space,
                       Kokkos::DefaultExecutionSpace::memory_space> DefaultDeviceType;
#endif
}
}

#endif
