cmake_minimum_required (VERSION 3.5)

project (compose CXX)
set (CMAKE_CXX_STANDARD 11)

function (prc var)
  message ("${var}: ${${var}}")
endfunction ()

find_package (MPI REQUIRED)

if (Kokkos_DIR)
  include (${Kokkos_DIR}/kokkos.cmake)
  set (Kokkos_INCLUDE ${Kokkos_DIR}/include)
else ()
  message (FATAL_ERROR "COMPOSE requires Kokkos_DIR")
endif ()

set (SOURCES
  cedr/cedr_caas.cpp
  cedr/cedr_local.cpp
  cedr/cedr_mpi.cpp
  cedr/cedr_qlt.cpp
  cedr/cedr_test.cpp
  cedr/cedr_test_1d_transport.cpp
  cedr/cedr_test_randomized.cpp
  cedr/cedr_util.cpp)

set (HEADERS
  cedr/cedr.hpp
  cedr/cedr_caas.hpp
  cedr/cedr_caas_inl.hpp
  cedr/cedr_cdr.hpp
  cedr/cedr_kokkos.hpp
  cedr/cedr_local.hpp
  cedr/cedr_local_inl.hpp
  cedr/cedr_mpi.hpp
  cedr/cedr_mpi_inl.hpp
  cedr/cedr_qlt.hpp
  cedr/cedr_qlt_inl.hpp
  cedr/cedr_test.hpp
  cedr/cedr_test_randomized.hpp
  cedr/cedr_util.hpp
  siqk/siqk.hpp
  siqk/siqk_defs.hpp
  siqk/siqk_geometry.hpp
  siqk/siqk_intersect.hpp
  siqk/siqk_quadrature.hpp
  siqk/siqk_search.hpp
  siqk/siqk_sqr.hpp)

if (NOT COMPOSE_TEST_MPIRUN)
  set (COMPOSE_TEST_MPIRUN mpirun)
endif ()
if (NOT COMPOSE_TEST_NRANK)
  set (COMPOSE_TEST_NRANK 8)
endif ()

set (COMPOSE_COMPILE_FLAGS "${MPI_COMPILE_FLAGS} ${KOKKOS_CXXFLAGS} ${CMAKE_CXX_FLAGS}")
set (COMPOSE_LINK_FLAGS "${MPI_LINK_FLAGS} ${KOKKOS_LDFLAGS}")
set (COMPOSE_INCLUDES "${Kokkos_INCLUDE}")
set (COMPOSE_LIBRARIES ${MPI_LIBRARIES} ${KOKKOS_LIBS})

prc(MPI_COMPILE_FLAGS)
prc(MPI_LINK_FLAGS)
prc(MPI_LIBRARIES)
add_library (${PROJECT_NAME} ${SOURCES})
set_target_properties (${PROJECT_NAME} PROPERTIES
  COMPILE_FLAGS ${COMPOSE_COMPILE_FLAGS}
  LINK_FLAGS ${COMPOSE_LINK_FLAGS})
target_include_directories (${PROJECT_NAME} PUBLIC cedr siqk)
target_include_directories (${PROJECT_NAME} PRIVATE siqk cedr)
target_include_directories (${PROJECT_NAME} PUBLIC ${COMPOSE_INCLUDES})
target_link_libraries (${PROJECT_NAME} ${COMPOSE_LIBRARIES})

install (TARGETS ${PROJECT_NAME} ARCHIVE DESTINATION lib)
install (FILES ${HEADERS} DESTINATION include/compose)

enable_testing ()
add_subdirectory(siqk)
add_subdirectory(cedr)
