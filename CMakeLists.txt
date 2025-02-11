cmake_minimum_required(VERSION 3.14)
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

if(POLICY CMP0069)
  # Actually enforce LTO when it is requested instead of silently ignoring it.
  cmake_policy(SET CMP0069 NEW)
endif()

if(POLICY CMP0110)
  # add_test() supports arbitrary characters in test names.
  cmake_policy(SET CMP0110 NEW)
endif()

# This path is used for FindXXX modules and common cmake function definitions.
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
enable_testing()

option(MCRL2_ENABLE_EXPERIMENTAL    "Enable compilation of experimental tools." OFF)
option(MCRL2_ENABLE_DEPRECATED      "Enable compilation of deprecated tools." OFF)
option(MCRL2_ENABLE_DEVELOPER       "Enable compilation of developer tools." OFF)
option(MCRL2_ENABLE_GUI_TOOLS       "Enable compilation of tools with a graphical user interface." ON)

option(MCRL2_ENABLE_DOCUMENTATION   "Enable generation of documentation." OFF)
option(MCRL2_ENABLE_DOC_MANUAL      "Enable generation of tool manual pages, which requires built tools." ON)
option(MCRL2_ENABLE_DOC_DOXYGEN     "Enable generation of code reference using Doxygen." ON)
option(MCRL2_ENABLE_DOC_PDFLATEX    "Enable generation of reference pdfs using pdflatex, which requires tikz and amsmath packages." ON)

option(MCRL2_ENABLE_BOOST_JSON_SUPPORT "Enable compilation of tools and libraries that require boost::json. This raises the minimal boost version requirement to 1.71.0." OFF)
option(MCRL2_ENABLE_SYLVAN          "Enable the Sylvan library required by the following symbolic tools: lpsreach, pbessolvesymbolic and ltsconvertsymbolic" ${UNIX})
option(MCRL2_ENABLE_BENCHMARKS      "Enable benchmarks. Build the 'benchmarks' target to generate the necessary files and tools. Run the benchmarks using ctest." OFF)

option(MCRL2_ENABLE_TESTS           "Enable generation of library and random test targets." OFF)
option(MCRL2_ENABLE_HEADER_TESTS    "Enable generation of headertest targets." OFF)
option(MCRL2_ENABLE_TOOL_TESTS      "Enable generation of tool tests." ON)
option(MCRL2_SKIP_LONG_TESTS        "Do not execute tests that take a long time to run." OFF)
option(MCRL2_TEST_JITTYC            "Also test the compiling rewriters in the library tests. This can be time consuming." OFF)

# Advanced options
option(MCRL2_ENABLE_STABLE           "Enable compilation of stable tools." ON)
option(MCRL2_ENABLE_ADDRESSSANITIZER "Enable additional memory run-time checks implemented by the AddressSanitizer." OFF)
option(MCRL2_ENABLE_MEMORYSANITIZER  "Enable additional memory run-time checks implemented by the MemorySanitizer." OFF)
option(MCRL2_ENABLE_THREADSANITIZER  "Enable additional data race run-time checks implemented by the ThreadSanitizer." OFF)
option(MCRL2_ENABLE_CODE_COVERAGE    "Enable generation of code coverage statistics." OFF)
option(MCRL2_ENABLE_STD_CHECKS       "Enable libstdc++ and libc++ specific assertions, e.g. safe iterator checks." OFF)
option(MCRL2_ENABLE_STD_CHECKS_BACKTRACE "Enables backtraces for the STD_CHECKS, only works on Ubuntu 24.04 and macOS" OFF)
option(MCRL2_ENABLE_SOUNDNESS_CHECKS "Enable extensive soundness checks." ON)
option(MCRL2_ENABLE_MSVC_CCACHE      "Enable MSVC specific flags to enable building with ccache" OFF)
option(MCRL2_ENABLE_JITTYC           "Enable the compiling rewriter (jittyc) functionality," ${UNIX})
option(MCRL2_ENABLE_MACHINENUMBERS   "Enable the usage of 64bit machine number digits to represent numeric sorts" ON)
option(MCRL2_ENABLE_MULTITHREADING   "Enable the usage of multiple threads. Disabling removes usage of synchronisation primitives" ON)
option(MCRL2_ENABLE_LINKER_LLD       "Enable LLVM lld as linker, which offers significantly linking speedup" OFF)

mark_as_advanced(
  MCRL2_ENABLE_ADDRESSSANITIZER
  MCRL2_ENABLE_MEMORYSANITIZER
  MCRL2_ENABLE_CODE_COVERAGE
  MCRL2_ENABLE_JITTYC
  MCRL2_ENABLE_LINKER_LLD
  MCRL2_ENABLE_MACHINENUMBERS
  MCRL2_ENABLE_MSVC_CCACHE
  MCRL2_ENABLE_MULTITHREADING
  MCRL2_ENABLE_SOUNDNESS_CHECKS 
  MCRL2_ENABLE_STABLE
  MCRL2_ENABLE_STD_CHECKS
  MCRL2_ENABLE_STD_CHECKS_BACKTRACE
  MCRL2_ENABLE_THREADSANITIZER
)

# A list of tools that is being built by the current configuration (populated by build/cmake/MCRL2AddTarget.cmake)
set_property(GLOBAL PROPERTY MCRL2_TOOLS)
set(MCRL2_QT_APPS "" CACHE INTERNAL "Internally keep track of Qt apps for the packaging procedure")
set(MCRL2_JITTYC_ARGUMENTS "" CACHE INTERNAL "Keeps track of compile flags used by the compiling rewriter")

if(MCRL2_ENABLE_BOOST_JSON_SUPPORT)
  set(MCRL2_MIN_BOOST_VERSION 1.75.0)
else()
  set(MCRL2_MIN_BOOST_VERSION 1.71.0)
endif()

if(MCRL2_ENABLE_GUI_TOOLS)
  find_package(OpenGL     QUIET REQUIRED)
  find_package(CVC3       QUIET)
  find_package(Qt6 6.2.4 COMPONENTS Core OpenGL OpenGLWidgets Widgets Xml Gui QUIET REQUIRED)
	
  # The directories for components are derived from Qt6_DIR, as such they can be marked as advanced.
  mark_as_advanced(Qt6_DIR Qt6Core_DIR Qt6OpenGL_DIR Qt6OpenGLWidgets_DIR Qt6Widgets_DIR Qt6Xml_DIR Qt6Gui_DIR)
endif()

find_package(Boost ${MCRL2_MIN_BOOST_VERSION} QUIET REQUIRED)

# We always require Threads even for single threaded application, but then threadsafe features are not used.
find_package(Threads)

include(ConfigurePlatform)
include(ConfigureCompiler)
include(MCRL2Version)
include(MCRL2AddTarget)

if(MCRL2_ENABLE_CODE_COVERAGE)
  include(CodeCoverage)
endif()

add_subdirectory(3rd-party/dparser)
if(MCRL2_ENABLE_GUI_TOOLS)
  add_subdirectory(3rd-party/tr)
endif()
if (MCRL2_ENABLE_SYLVAN)
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

# Install the example directory.
install(DIRECTORY examples DESTINATION ${MCRL2_RESOURCE_PATH} COMPONENT Examples)

# Install third party libraries if they are shared libraries
install(TARGETS dparser
  COMPONENT "Libraries"
  LIBRARY DESTINATION ${MCRL2_LIBRARY_PATH}
  ARCHIVE DESTINATION ${MCRL2_ARCHIVE_PATH}
  FRAMEWORK DESTINATION ${MCRL2_LIBRARY_PATH})

if(MCRL2_ENABLE_GUI_TOOLS)
  install(TARGETS tr
    COMPONENT "Libraries"
    LIBRARY DESTINATION ${MCRL2_LIBRARY_PATH}
    ARCHIVE DESTINATION ${MCRL2_ARCHIVE_PATH}
    FRAMEWORK DESTINATION ${MCRL2_LIBRARY_PATH})
endif()

if (MCRL2_ENABLE_SYLVAN)
  install(TARGETS sylvan
    COMPONENT "Libraries"
    LIBRARY DESTINATION ${MCRL2_LIBRARY_PATH}
    ARCHIVE DESTINATION ${MCRL2_ARCHIVE_PATH}
    FRAMEWORK DESTINATION ${MCRL2_LIBRARY_PATH})
endif()

# Add the packaging script, done after all the configuration has finished.
add_subdirectory(cmake/packaging)

# Prints dependency versions to the cmake configure output.
include(PrintBuildInfo)
