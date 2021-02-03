cmake_minimum_required(VERSION 3.8)
project(mCRL2)

if(POLICY CMP0054)
  cmake_policy(SET CMP0054 NEW)
endif()

if(POLICY CMP0065)
  # Do not export symbols from executables unless ENABLE_EXPORTS is explicitly set.
  cmake_policy(SET CMP0065 NEW)
endif()

if(POLICY CMP0071)
  # The AUTOMOC and AUTOUIC properties should not apply to GENERATED files.
  cmake_policy(SET CMP0071 NEW)
endif()

if(POLICY CMP0072)
  # Choose the GLVND (GL Vendor Neutral Dispatch library) over OpenGL legacy if both are available.
  cmake_policy(SET CMP0072 NEW)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/build/cmake)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
enable_testing()

if(APPLE)
  # Lower versions of Qt do not support native menu's correctly on Apple.
  # To use a lower version replace set window->menuBar()->setNativeMenuBar(true)
  # by window->menuBar()->setNativeMenuBar(false) in line 138 of libraries/gui/include/mcrl2/gui/qt_tool.h.
  set(MCRL2_MIN_QT_VERSION 5.10.0)
else()
  set(MCRL2_MIN_QT_VERSION 5.9.0)
endif()

if (UNIX)
  # On Ubuntu bionic the available libboost version is the following:
  set(MCRL2_MIN_BOOST_VERSION 1.65.1)
else()
  set(MCRL2_MIN_BOOST_VERSION 1.66.0)
endif()

option(MCRL2_ENABLE_GUI_TOOLS       "Enable compilation of tools with a graphical user interface." ON)
option(MCRL2_ENABLE_EXPERIMENTAL    "Enable compilation of experimental tools." OFF)
option(MCRL2_ENABLE_DEPRECATED      "Enable compilation of deprecated tools." OFF)
option(MCRL2_ENABLE_DEVELOPER       "Enable compilation of developer tools." OFF)
option(MCRL2_ENABLE_HEADER_TESTS    "Enable generation of headertest targets." OFF)
option(MCRL2_ENABLE_TESTS           "Enable generation of library, tool and random test targets." OFF)
option(MCRL2_ENABLE_DOCUMENTATION   "Enable generation of documentation." OFF)
option(MCRL2_ENABLE_BENCHMARKS      "Enable benchmarks. Build the 'benchmarks' target to generate the necessary files and tools. Run the benchmarks using ctest." OFF)
option(MCRL2_EXTRA_TOOL_TESTS       "Enable testing of tools on more mCRL2 specifications." OFF)
option(MCRL2_TEST_JITTYC            "Also test the compiling rewriters in the library tests. This can be time consuming." OFF)
set(MCRL2_QT_APPS "" CACHE INTERNAL "Internally keep track of Qt apps for the packaging procedure")

option(MCRL2_ENABLE_ADDRESSSANITIZER "Enable additional memory run-time checks implemented by the AddressSanitizer." OFF)
option(MCRL2_ENABLE_CODE_COVERAGE   "Enable generation of code coverage statistics." OFF)
option(MCRL2_ENABLE_DEBUG_SOUNDNESS_CHECKS "Enable extensive soundness check in the Debug build type." ON)
option(MCRL2_ENABLE_STABLE          "Enable compilation of stable tools." ON)
option(MCRL2_SKIP_LONG_TESTS        "Do not execute tests that take a long time to run." OFF)

mark_as_advanced(
  MCRL2_ENABLE_ADDRESSSANITIZER
  MCRL2_ENABLE_CODE_COVERAGE
  MCRL2_ENABLE_DEBUG_SOUNDNESS_CHECKS 
  MCRL2_ENABLE_STABLE
)

if(MCRL2_ENABLE_GUI_TOOLS)
  find_package(OpenGL     QUIET REQUIRED)
  find_package(TR         QUIET REQUIRED)
  find_package(CVC3       QUIET)
  find_package(Qt5 ${MCRL2_MIN_QT_VERSION} COMPONENTS Core OpenGL Widgets Xml Gui QUIET REQUIRED)
	
  # The directories for components are derived from Qt5_DIR, as such they can be marked as advanced.
  mark_as_advanced(Qt5_DIR Qt5Core_DIR Qt5OpenGL_DIR Qt5Widgets_DIR Qt5Xml_DIR Qt5Gui_DIR)
endif()

find_package(Boost ${MCRL2_MIN_BOOST_VERSION} QUIET REQUIRED)

include(ConfigurePlatform)
include(ConfigureCompiler)
include(MCRL2Version)
include(AddMCRL2Binary)

if(MCRL2_ENABLE_CODE_COVERAGE)
  include(CodeCoverage)
endif()

add_subdirectory(3rd-party/dparser)

if (UNIX)
  add_subdirectory(3rd-party/sylvan/src)
endif()

if (MCRL2_ENABLE_BENCHMARKS)
  add_subdirectory(benchmarks)
endif()

add_subdirectory(libraries)
add_subdirectory(tools)

if (MCRL2_ENABLE_DOCUMENTATION)
  add_subdirectory(doc)
endif()

if (MCRL2_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

install(DIRECTORY examples DESTINATION ${MCRL2_RESOURCE_PATH} COMPONENT Examples)

# Add the packaging as a subdirectory to assure that CMake only performs install
# commands from this script _after_ all install commands from other subdirectories
# are completed.
add_subdirectory(build/packaging)

include(PrintBuildInfo)
