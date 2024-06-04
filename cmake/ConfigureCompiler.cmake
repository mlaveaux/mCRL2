# Perform compiler-specific compiler configuration
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 19.20)
    message(FATAL_ERROR "MSVC version 19.20 (aka MSVC 2019 v16.0) is required.")
  endif()
  include(ConfigureMSVC)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

  # We do not test on GCC9 directly, but it used on the mastodont.
  if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 11.0)
    message(FATAL_ERROR "GCC version must be at least 11.0.")
  endif()
  include(ConfigureUNIX)

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

  if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 14.0)
    message(FATAL_ERROR "Clang version must be at least 14.0.")
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
if(MCRL2_ENABLE_OXIDD)
  add_compile_definitions(MCRL2_ENABLE_OXIDD)
endif()


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
endif(MCRL2_SKIP_LONG_TESTS)

if(APPLE)
  # Silence useless OpenGL deprecration warnings on macOS. Some GUI tools use outdated OpenGL and this will 
  # only be replaced when it is removed.
  add_compile_definitions(GL_SILENCE_DEPRECATION)
endif()

# Enable C++17 for all targets.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED true)
