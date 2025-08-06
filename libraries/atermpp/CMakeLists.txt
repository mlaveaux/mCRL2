mcrl2_add_library(mcrl2_atermpp
  SOURCES
    source/aterm_implementation.cpp
    source/aterm_io_binary.cpp
    source/aterm_io_text.cpp
    source/function_symbol.cpp
    source/function_symbol_pool.cpp
  DEPENDS
    mcrl2_utilities
    Threads::Threads
)

if (MCRL2_ENABLE_MODULES)
  target_sources(mcrl2_atermpp
    PUBLIC
      FILE_SET CXX_MODULES FILES
        include/mcrl2/atermpp/aterm.cxx
        include/mcrl2/atermpp/aterm_core.cxx
        include/mcrl2/atermpp/aterm_list.cxx
        include/mcrl2/atermpp/concepts.cxx
        include/mcrl2/atermpp/type_traits.cxx
        include/mcrl2/atermpp/function_symbol.cxx
        include/mcrl2/atermpp/module.cxx
        include/mcrl2/atermpp/detail/aterm_configuration.cxx
        include/mcrl2/atermpp/detail/aterm_container.cxx
        include/mcrl2/atermpp/detail/thread_aterm_pool.cxx
        include/mcrl2/atermpp/detail/aterm_pool_storage.cxx
        include/mcrl2/atermpp/detail/aterm_pool.cxx
        include/mcrl2/atermpp/detail/global_aterm_pool.cxx
        include/mcrl2/atermpp/detail/aterm_data.cxx
        include/mcrl2/atermpp/detail/aterm_core_data.cxx
        include/mcrl2/atermpp/detail/aterm_list_data.cxx
        include/mcrl2/atermpp/detail/aterm_list_iterator.cxx
        include/mcrl2/atermpp/detail/aterm_hash.cxx
        include/mcrl2/atermpp/detail/aterm_int_data.cxx
        include/mcrl2/atermpp/detail/index_traits.cxx
        include/mcrl2/atermpp/detail/function_symbol_data.cxx
        include/mcrl2/atermpp/detail/function_symbol_hash.cxx
        include/mcrl2/atermpp/detail/function_symbol_pool.cxx
        include/mcrl2/atermpp/detail/aterm_appl_iterator.cxx
        include/mcrl2/atermpp/standard_containers/unordered_map.cxx
  )
endif()


if (${MCRL2_ENABLE_BENCHMARKS})
  add_subdirectory(benchmark/)
endif()

add_subdirectory(example)
