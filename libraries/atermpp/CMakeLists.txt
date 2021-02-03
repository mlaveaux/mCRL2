add_mcrl2_library(atermpp
  INSTALL_HEADERS TRUE
  SOURCES
    aterm_implementation.cpp
    aterm_io_binary.cpp
    aterm_io_text.cpp
    function_symbol.cpp
    function_symbol_pool.cpp
  DEPENDS
    mcrl2_utilities
)

if (${MCRL2_ENABLE_BENCHMARKS})
  add_subdirectory(benchmark/)
endif()

add_subdirectory(example)
