set(SPHINX_BUILD_CMD "sphinx-build")
set(SPHINX_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR})
set(SPHINX_BUILD_ROOT_DIR "${CMAKE_BINARY_DIR}/sphinx")
set(SPHINX_BUILD_SRC_DIR "${SPHINX_BUILD_ROOT_DIR}/source")
set(SPHINX_BUILD_OUT_DIR "${SPHINX_BUILD_ROOT_DIR}/build")
set(SPHINX_BUILD_TEMP_DIR "${SPHINX_BUILD_ROOT_DIR}/temp")

set(SPHINX_PARALLEL "-jauto")
if(CMAKE_BUILD_PARALLEL_LEVEL)
  set(SPHINX_PARALLEL "-j${CMAKE_BUILD_PARALLEL_LEVEL}")
endif()

add_custom_target(cleandoc
  COMMENT "Delete all cached documentation"
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${SPHINX_BUILD_SRC_DIR}
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${SPHINX_BUILD_OUT_DIR}
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${SPHINX_BUILD_TEMP_DIR}
  VERBATIM
)

if(MCRL2_ENABLE_DOC_DOXYGEN)
  # Replace variables @@ with the current values
  set(DOXYGEN_BASE ${CMAKE_BINARY_DIR}/sphinx/source/_doxygen/)
  set(DOXYGEN_OUTPUT_PATH ${DOXYGEN_BASE}/output/)
  set(DOXYGEN_TAG_PATH ${DOXYGEN_BASE}/mcrl2.tag)
  set(DOXYGEN_INPUT_PATH ${CMAKE_SOURCE_DIR}/libraries/)

  configure_file(${CMAKE_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_BINARY_DIR}/doc/Doxyfile @ONLY)

  # Configure two doxygen commands (one will have the cleandoc as dependency)
  foreach(NAME "doxygen" "doxygen_clean")
    add_custom_target(${NAME}
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/libraries/
      COMMENT "Extracting source code documentation using Doxygen"
      # ensure the build directory exists, because Doxygen does not create it
      # automatically
      COMMAND ${CMAKE_COMMAND} -E make_directory ${DOXYGEN_OUTPUT_PATH}
      COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/doc/Doxyfile
      VERBATIM
      USES_TERMINAL)
  endforeach()
endif()

# Get the list of tools that are built by the current configuration
get_property(MCRL2_TOOLS GLOBAL PROPERTY MCRL2_TOOLS)

function(add_build_doc_target NAME DESC)
  set(SPHINX_TAGS "")
  if(MCRL2_ENABLE_DOC_DOXYGEN)
    list(APPEND SPHINX_TAGS -t build_doxygen)
  endif()
  if(MCRL2_ENABLE_DOC_PDFLATEX)
    list(APPEND SPHINX_TAGS -t build_pdflatex)
  endif()
  if(MCRL2_ENABLE_DOC_MANUAL)
    list(APPEND SPHINX_TAGS -t build_manual)
  endif()
  
  # Replace the cmake variables by their values (we do not use the environment since we want to be able to call sphinx_build from the command line)
  configure_file(${SPHINX_SOURCE_DIR}/conf.in.py ${SPHINX_BUILD_ROOT_DIR}/conf.py @ONLY)

  add_custom_target(${NAME}
    COMMENT ${DESC}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    # copy the Sphinx source to the build directory, because some of the
    # extensions modify the source and we don't want to pollute the source dir
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${SPHINX_SOURCE_DIR} ${SPHINX_BUILD_SRC_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy ${SPHINX_BUILD_ROOT_DIR}/conf.py ${SPHINX_BUILD_SRC_DIR}/conf.py
    # run the Sphinx build command
    COMMAND ${SPHINX_BUILD_CMD}
      -b html 
      -d "${SPHINX_BUILD_TEMP_DIR}/doctrees"
      -T
#     -W --keep-going
      ${SPHINX_TAGS}
      ${SPHINX_PARALLEL}
      ${SPHINX_BUILD_SRC_DIR} 
      ${SPHINX_BUILD_OUT_DIR}      
    VERBATIM
    USES_TERMINAL
  )
endfunction()

add_build_doc_target(fastdoc "Generating mCRL2 documentation (does not clean up old files)")
add_build_doc_target(doc "Generating mCRL2 documentation")

if (MCRL2_ENABLE_DOC_DOXYGEN)
  add_dependencies(fastdoc doxygen)

  # Ensure that cleandoc is called before doxygen_clean
  add_dependencies(doc doxygen_clean)
  add_dependencies(doxygen_clean cleandoc)
else()
  add_dependencies(doc cleandoc)
endif()

if(MCRL2_ENABLE_DOCUMENTATION_MANUAL)
  if(NOT MCRL2_ENABLE_EXPERIMENTAL or NOT MCRL2_ENABLE_STABLE)
    message(FATAL_ERROR "The option MCRL2_ENABLE_DOCUMENTATION_MANUAL requires both the stable and experimental tools to be built")
  endif()
endif()