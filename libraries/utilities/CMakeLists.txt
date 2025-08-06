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
        include/mcrl2/utilities/block_allocator.cxx
        include/mcrl2/utilities/cache_metric.cxx
        include/mcrl2/utilities/configuration.cxx
        include/mcrl2/utilities/detail/atomic_wrapper.cxx
        include/mcrl2/utilities/detail/bucket_list.cxx
        include/mcrl2/utilities/detail/free_list.cxx
        include/mcrl2/utilities/exception.cxx
        include/mcrl2/utilities/execution_timer.cxx
        include/mcrl2/utilities/hash_utility.cxx
        include/mcrl2/utilities/hashtable.cxx
        include/mcrl2/utilities/indexed_set.cxx
        include/mcrl2/utilities/logger.cxx
        include/mcrl2/utilities/memory_pool.cxx
        include/mcrl2/utilities/module.cxx
        include/mcrl2/utilities/mutex.cxx  
        include/mcrl2/utilities/noncopyable.cxx   
        include/mcrl2/utilities/power_of_two.cxx
        include/mcrl2/utilities/shared_mutex.cxx  
        include/mcrl2/utilities/shared_reference.cxx
        include/mcrl2/utilities/stack_array.cxx
        include/mcrl2/utilities/tagged_pointer.cxx
        include/mcrl2/utilities/text_utility.cxx
        include/mcrl2/utilities/type_traits.cxx     
        include/mcrl2/utilities/unordered_map.cxx
        include/mcrl2/utilities/unordered_set.cxx
  )
endif()

add_subdirectory(example)
