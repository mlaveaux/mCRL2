

mcrl2_add_library(mcrl2_symbolic
  SOURCES
    source/ldd_stream.cpp
  DEPENDS
    mcrl2_data
    Boost::boost
)

if(MCRL2_ENABLE_OXIDD)
  target_link_libraries(mcrl2_symbolic PUBLIC oxidd)
endif()

if(MCRL2_ENABLE_SYLVAN)
  # Only depend on Sylvan when we have enabled Sylvan compilation.
  target_link_libraries(mcrl2_symbolic PUBLIC sylvan)
endif()
