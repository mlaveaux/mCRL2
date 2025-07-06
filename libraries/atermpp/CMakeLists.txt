mcrl2_add_library(mcrl2_atermpp
  SOURCES
    source/aterm_io_binary.cpp
    source/aterm_io_text.cpp
  DEPENDS
    mcrl2_utilities
    Threads::Threads
)

target_link_libraries(mcrl2_atermpp PUBLIC mcrl3_ffi)

if (${MCRL2_ENABLE_BENCHMARKS})
  add_subdirectory(benchmark/)
endif()

add_subdirectory(example)
