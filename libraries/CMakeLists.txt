###
### The mCRL2 libraries
###

set(MCRL2_LIBRARIES
  atermpp
  bes
  core
  data
  gui
  lps
  lts
  modal_formula
  pbes
  pg
  process
  smt
  utilities
)

foreach(LIB ${MCRL2_LIBRARIES})
  add_subdirectory(${LIB})
endforeach()

# Add custom directory if it exists.
IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/custom)
  add_subdirectory(custom)
endif()
