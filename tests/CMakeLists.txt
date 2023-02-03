include(GenerateToolTests)
find_package(PythonModules COMPONENTS yaml psutil REQUIRED)

generate_tool_tests(${CMAKE_SOURCE_DIR}/examples/academic/abp/abp.mcrl2 "")

if(MCRL2_EXTRA_TOOL_TESTS)
  set(MCRL2_DIR ${CMAKE_SOURCE_DIR}/examples/language)
  generate_tool_tests(${MCRL2_DIR}/list.mcrl2 "lpsconfcheck")
  generate_tool_tests(${MCRL2_DIR}/funccomp.mcrl2 "lpsparunfold")
  # lpsrealelm cannot deal with "P / 2" in the following
  generate_tool_tests(${MCRL2_DIR}/divide2_10.mcrl2 "lpsrealelm")
  # Disable lts2*, because LTS contains a "Terminate" actions, which
  # is not declared in the mCRL2 spec
  generate_tool_tests(${MCRL2_DIR}/tau.mcrl2 "lts2lps;lts2pbes")
  generate_tool_tests(${MCRL2_DIR}/par.mcrl2 "lts2lps;lts2pbes")
  # Disable lpsparunfold when LPS has no process parameters
  generate_tool_tests(${MCRL2_DIR}/struct.mcrl2 "lpsparunfold")
  generate_tool_tests(${MCRL2_DIR}/forall.mcrl2 "lpsparunfold")
  generate_tool_tests(${MCRL2_DIR}/exists.mcrl2 "lpsparunfold")
  generate_tool_tests(${MCRL2_DIR}/small3.mcrl2 "lpsparunfold")
  generate_tool_tests(${MCRL2_DIR}/small2.mcrl2 "lpsparunfold")
  generate_tool_tests(${MCRL2_DIR}/upcast.mcrl2 "lpsparunfold")
  generate_tool_tests(${MCRL2_DIR}/delta0.mcrl2 "lpsparunfold")
  generate_tool_tests(${MCRL2_DIR}/lambda.mcrl2 "lpsparunfold")
  generate_tool_tests(${MCRL2_DIR}/delta.mcrl2 "lpsparunfold")
  generate_tool_tests(${MCRL2_DIR}/time.mcrl2 "lpsparunfold;lpsrealelm")
  # LPS and PBES instantiation time out on most platforms. All tests
  # that depend on the LTS or BES formats are also disabled.
  generate_tool_tests(${MCRL2_DIR}/sets_bags.mcrl2 "lpsbinary;lts2lps;lts2pbes;lps2lts;ltsinfo;ltsconvert;ltscompare;pbes2bool;pbes2bes;pbesinst;besconvert;bespp;besinfo;txt2bes;bessolve;pbespgsolve")
  generate_tool_tests(${MCRL2_DIR}/rational.mcrl2 "pbes2bool;pbesrewr;lts2lps;besinfo;bespp;lpsbinary;pbes2bes;besconvert;bessolve;pbesinst;lts2pbes;pbespgsolve;txt2bes;lpsrealelm;lps2lts;ltsinfo;ltsconvert;ltscompare")
endif()

# The workspace is used for temporary files (it is the directory where the test script is executed from).
set(WORKSPACE ${CMAKE_CURRENT_BINARY_DIR}/workspace)
file(MAKE_DIRECTORY ${WORKSPACE})

#
# Python test drivers are assumed to have the following behaviour random_testing.py --names
# prints all test names (on separate lines) that the test script can run.
# For each <name> a target is created named random_<name> that runs 100 consecutive
# tests with --repetitions 100 choosing the tests with the --pattern <name> option.
#
foreach(type "random" "regression")
  set(testrunner ${CMAKE_CURRENT_SOURCE_DIR}/${type}/${type}_testing.py)
  execute_process(COMMAND ${PYTHON_EXECUTABLE} ${testrunner} --names
    WORKING_DIRECTORY ${WORKSPACE}
    OUTPUT_VARIABLE available_tests
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE testnames_error)
  if(testnames_error)
    message(SEND_ERROR "${testrunner} did not succesfully list its test targets. Please contact an mCRL2 developer to fix this.")
  endif()

  string(REPLACE "\n" ";" available_tests "${available_tests}")
  string(REPLACE "\r" "" available_tests "${available_tests}")

  foreach(testname ${available_tests})
    set(REPETITIONS "")
    if("${type}" STREQUAL "random")
      set(REPETITIONS "-r 100")
    endif()

    set(targetname "${type}_${testname}")

    add_test(NAME ${targetname}
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/workspace
             COMMAND ${PYTHON_EXECUTABLE} ${testrunner} ${REPETITIONS} --toolpath ${MCRL2_STAGE_PATH}/${MCRL2_RUNTIME_PATH} --pattern ${testname})
    set_tests_properties(${targetname} PROPERTIES LABELS "${type}test")
  endforeach()
endforeach()
