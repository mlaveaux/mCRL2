# mcrl2server uses Boost.JSON, which is only available from Boost 1.75 onwards.
if(Boost_VERSION VERSION_LESS 1.75.0)
  message(STATUS "Not building mcrl2server: Boost.JSON requires Boost 1.75 or newer (found ${Boost_VERSION}).")
  return()
endif()

mcrl2_add_tool(mcrl2server
  SOURCES
    mcrl2server.cpp
  DEPENDS
    mcrl2_lps
    Boost::boost
    Threads::Threads
)
