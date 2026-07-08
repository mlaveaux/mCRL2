# Authors: mCRL2 developers
# Copyright: see the accompanying file COPYING or copy at
# https://github.com/mCRL2org/mCRL2/blob/master/COPYING
#
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

# Compile (and run) the C++ example snippets that are embedded in the developer
# manual via `literalinclude`. Building them as a "doctest" suite ensures the
# documented code keeps compiling and behaving as described.
file(GLOB DOCTEST_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

foreach(source ${DOCTEST_SOURCES})
  get_filename_component(base ${source} NAME_WE)
  set(testname doctest_core_${base})

  add_executable(${testname} ${source})
  set_target_properties(${testname} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/doctest)
  target_link_libraries(${testname} mcrl2_core)

  add_test(NAME ${testname}
    COMMAND "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/doctest/${testname}"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/Testing/")
  set_tests_properties(${testname} PROPERTIES LABELS doctest)
endforeach()
