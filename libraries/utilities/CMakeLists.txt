# Install a header file containing the version number
configure_file(include/mcrl2/utilities/toolset_version_const.h.in mcrl2/utilities/toolset_version_const.h)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/mcrl2/utilities/toolset_version_const.h
        DESTINATION ${MCRL2_INCLUDE_PATH}/mcrl2/utilities
        COMPONENT Headers)

mcrl2_add_library(mcrl2_utilities
  SOURCES
    source/bitstream.cpp
    source/cache_metric.cpp
    source/command_line_interface.cpp
    source/logger.cpp
    source/text_utility.cpp
    source/toolset_version.cpp
  INCLUDE_DIRS
    ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS
    Boost::boost
    Threads::Threads
)

if (MCRL2_ENABLE_MODULES)
  target_sources(mcrl2_utilities
    PUBLIC
      FILE_SET CXX_MODULES FILES
        include/mcrl2/utilities/utilities.cxx
        include/mcrl2/utilities/exception.cxx
        include/mcrl2/utilities/execution_timer.cxx
  )
endif()

add_subdirectory(example)
