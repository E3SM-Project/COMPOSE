add_executable (siqk_test siqk_test.cpp)
set_target_properties (siqk_test PROPERTIES
  COMPILE_FLAGS ${COMPOSE_COMPILE_FLAGS}
  LINK_FLAGS ${COMPOSE_LINK_FLAGS})
target_include_directories (siqk_test PRIVATE ${COMPOSE_INCLUDES})
target_link_libraries (siqk_test ${COMPOSE_LIBRARIES})

configure_file (siqk_runtests.py siqk_runtests.py)

find_package(Python COMPONENTS Interpreter)

if (Python_FOUND)
add_test (NAME siqk-test-area
  COMMAND ${Python_EXECUTABLE} siqk_runtests.py $<TARGET_FILE:siqk_test> 0)
add_test (NAME siqk-test-cube
  COMMAND ${Python_EXECUTABLE} siqk_runtests.py $<TARGET_FILE:siqk_test> 1)
else()
message(STATUS "python not found; disabling siqk tests.")
endif()
