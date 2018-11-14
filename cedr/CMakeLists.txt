add_executable (cedr_test cedr_test.cpp)
set_target_properties (cedr_test PROPERTIES
  COMPILE_FLAGS ${COMPOSE_COMPILE_FLAGS}
  LINK_FLAGS ${COMPOSE_LINK_FLAGS})

target_include_directories (cedr_test PRIVATE ${COMPOSE_INCLUDES})
target_link_libraries (cedr_test ${PROJECT_NAME} ${COMPOSE_LIBRARIES})

add_test (NAME cedr-test-unit
  COMMAND $<TARGET_FILE:cedr_test> -t)
add_test (NAME cedr-test-unit-mpi
  COMMAND ${COMPOSE_TEST_MPIRUN} ${COMPOSE_TEST_MPIFLAGS} -np ${COMPOSE_TEST_NRANK}
  $<TARGET_FILE:cedr_test> -t --proc-random -nc 111 -nt 11)
add_test (NAME cedr-test-t1d
  COMMAND $<TARGET_FILE:cedr_test> -t -t1d -nc 111)
