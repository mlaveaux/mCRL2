# Perform compiler-specific compiler configuration
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 19.31)
    message(FATAL_ERROR "MSVC version 19.31 (aka MSVC 2022 v17.1) is required.")
  endif()
  include(ConfigureMSVC)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

  # We do not test on GCC9 directly, but it used on the mastodont.
  if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 11.0)
    message(FATAL_ERROR "GCC version must be at least 11.0.")
  endif()
  include(ConfigureUNIX)

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

  if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 16.0)
    message(FATAL_ERROR "Clang version must be at least 16.0.")
  endif()
  set(MCRL2_IS_CLANG 1)
  include(ConfigureUNIX)

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  # The following should actually be version 10.0.1, but then VERSION_LESS does not handle AppleClang 11 correctly.
  if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 14.0)
    message(FATAL_ERROR "AppleClang version must be at least 14.0")
  endif()
  include(ConfigureUNIX)
else()
  message(FATAL_ERROR "Unsupported compiler setup (C: ${CMAKE_C_COMPILER_ID} / C++: ${CMAKE_CXX_COMPILER_ID}).")
endif()

# Sets the default build type if none was provided.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Sets the build type, defaults to Release." FORCE)
endif()

if(CMAKE_BUILD_TYPE)
  # Enables a selection of build types instead of typing it.
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo") 
endif()

# Add the definition to disable soundness checks when the configuration is set to OFF.
if(NOT ${MCRL2_ENABLE_SOUNDNESS_CHECKS})
  add_compile_definitions(MCRL2_NO_SOUNDNESS_CHECKS)
endif()

# Add compiler definitions
if(MCRL2_ENABLE_SYLVAN)
  add_compile_definitions(MCRL2_ENABLE_SYLVAN)
endif()

if(MCRL2_ENABLE_JITTYC)
  add_compile_definitions(MCRL2_ENABLE_JITTYC)

  if(MCRL2_TEST_JITTYC)
    add_compile_definitions(MCRL2_TEST_JITTYC)
  endif()
endif()

if(MCRL2_ENABLE_MACHINENUMBERS)
  add_compile_definitions(MCRL2_ENABLE_MACHINENUMBERS)
endif()

if(MCRL2_ENABLE_MULTITHREADING)
  add_compile_definitions(MCRL2_ENABLE_MULTITHREADING)
endif()

if(MCRL2_SKIP_LONG_TESTS)
  add_compile_definitions(MCRL2_SKIP_LONG_TESTS)
endif()

if(APPLE)
  # Silence useless OpenGL deprecration warnings on macOS. Some GUI tools use outdated OpenGL and this will 
  # only be replaced when it is removed.
  add_compile_definitions(GL_SILENCE_DEPRECATION)
endif()

if(MCRL2_ENABLE_MODULES)
  add_compile_definitions(MCRL2_ENABLE_MODULES)
  
  # Only add the export keyword when modules are enabled, otherwise leave blank.
  add_compile_definitions(MCRL2_EXPORT=export)
endif()

# Enable the new JFG branching bisimulation algorithm.
add_compile_definitions(BRANCH_BIS_EXPERIMENT_JFG)

# Enable all macros defined in the code used for debugging purposes.
add_debug_compile_definitions(MCRL2_PBES_STATEGRAPH_CHECK_GUARDS)


# Enable all macros defined in the code used for debugging purposes.
add_debug_compile_definitions(MCRL2_PBES_STATEGRAPH_CHECK_GUARDS)
add_debug_compile_definitions(MCRL2_LPS_PARELM_DEBUG)
add_debug_compile_definitions(MCRL2_ABSINTHE_CHECK_EXPRESSIONS)
add_debug_compile_definitions(PARANOID_CHECK)
add_compile_definitions(MCRL2_EXTENDED_TESTS)

# These are defines that can be enabled for additional debug printing
#add_compile_definitions(MCRL2_DEBUG_EXPRESSION_BUILDER)
#add_compile_definitions(MCRL2_PBES_EXPRESSION_BUILDER_DEBUG)
#add_compile_definitions(MCRL2_PFNF_VISITOR_DEBUG)
#add_compile_definitions(MCRL2_SMALL_PROGRESS_MEASURES_DEBUG)
#add_compile_definitions(MCRL2_LOG_ENUMERATOR)
#add_compile_definitions(MCRL2_ENUMERATOR_COUNT_REWRITE_CALLS)
#add_compile_definitions(MCRL2_PRINT_REWRITE_STEPS)
#add_compile_definitions(MCRL2_EQUAL_MULTI_ACTIONS_DEBUG)
#add_compile_definitions(MCRL2_LOG_LPS_LINEARISE_STATISTICS)
#add_compile_definitions(MCRL2_BES2PGSOLVER_PRINT_VARIABLE_NAMES)

# TODO: These belong the the dnj algorithm.
#add_compile_definitions(TEST_WORK_COUNTER_NAMES)
#add_compile_definitions(USE_SIMPLE_LIST)
#add_compile_definitions(USE_POOL_ALLOCATOR)
#add_compile_definitions(INIT_WITHOUT_BLC_SETS)

# TODO: Is this still needed?
#add_compile_definitions(MCRL2_PGSOLVER_ENABLED)

# Enable C++20 for all targets.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)