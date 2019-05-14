if(NOT WIN32)
  set(COMPILING_REWRITER_SRC detail/rewrite/jittyc.cpp)
  set(COMPILING_REWRITER_DEPS dl)
  include(CompilingRewriter.cmake)
endif()

add_mcrl2_library(data
  INSTALL_HEADERS TRUE
  NOHEADERTEST
    mcrl2/data/detail/rewrite/jittycpreamble.h
  SOURCES
    data.cpp
    data_io.cpp
    data_specification.cpp
    typecheck.cpp
    detail/prover/smt_lib_solver.cpp
    detail/rewrite/innermost.cpp
    detail/rewrite/jitty.cpp
    detail/rewrite/rewrite.cpp
    detail/rewrite/strategy.cpp
    detail/match/naive_matcher.cpp
    ${COMPILING_REWRITER_SRC}
  DEPENDS
    mcrl2_core
    mcrl2_utilities
    ${COMPILING_REWRITER_DEPS}
)

add_subdirectory(example)
