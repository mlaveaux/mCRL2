# mcrl2server and mcrl2mbt use Boost.JSON, which is only available from Boost 1.75 onwards.
if(Boost_VERSION VERSION_LESS 1.75.0)
  message(STATUS "Not building mcrl2server/mcrl2mbt: Boost.JSON requires Boost 1.75 or newer (found ${Boost_VERSION}).")
  return()
endif()

mcrl2_add_tool(mcrl2mbt
  SOURCES
    mcrl2mbt.cpp
    io_classifier.cpp
    mbt_protocol.cpp
    mbt_client.cpp
    mbt_session.cpp
    state_set.cpp
  DEPENDS
    mcrl2_lps
    Boost::boost
    Threads::Threads
)

if(MCRL2_ENABLE_TESTS)
  macro(mcrl2mbt_add_test NAME)
    add_executable(${NAME} test/${NAME}.cpp)
    target_link_libraries(${NAME} mcrl2_lps Boost::boost Threads::Threads)
    target_include_directories(${NAME} PRIVATE ".")
    add_test(NAME ${NAME} COMMAND ${NAME})
  endmacro()

  mcrl2mbt_add_test(io_classifier_test)
  mcrl2mbt_add_test(mbt_protocol_test)
  mcrl2mbt_add_test(state_set_test)
  mcrl2mbt_add_test(mock_adapter_test)
endif()
