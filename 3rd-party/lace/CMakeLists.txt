cmake_minimum_required(VERSION 3.23...3.31)
cmake_policy(SET CMP0092 NEW) # do not automatically add /W3 for MSVC

project(lace
    VERSION 1.6.3
    DESCRIPTION "Lace, a work-stealing framework for multi-core fork-join parallelism"
    HOMEPAGE_URL "https://github.com/trolando/lace"
    LANGUAGES C
)

# ── Build defaults ───────────────────────────────────────────────────────────

if(lace_IS_TOP_LEVEL AND NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

if(lace_IS_TOP_LEVEL)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# ── Options ──────────────────────────────────────────────────────────────────

option(BUILD_SHARED_LIBS "Build shared libraries instead of static"                  OFF)
option(LACE_ENABLE_PIC   "Build Lace with -fPIC (e.g. for linking into shared libs)" OFF)
option(LACE_NATIVE_OPT   "Enable -march=native for compiler optimizations"           ON)
option(LACE_BACKOFF      "Let idle threads sleep if there is no work"                ON)
option(LACE_USE_HWLOC    "Let Lace pin threads/memory using libhwloc"                OFF)
option(LACE_PIE_TIMES    "Let Lace record Pie times"                                 OFF)
option(LACE_COUNT_TASKS  "Let Lace record the number of tasks"                       OFF)
option(LACE_COUNT_STEALS "Let Lace count #steals and #leaps"                         OFF)
option(LACE_COUNT_SPLITS "Let Lace count #splits"                                    OFF)

if(lace_IS_TOP_LEVEL)
    option(LACE_BUILD_TESTS      "Build and run Lace tests"                          ON)
    option(LACE_BUILD_BENCHMARKS "Build Lace benchmark programs"                     ON)
    option(LACE_SANITIZE_ADDRESS "Enable AddressSanitizer"                           OFF)
    option(LACE_SANITIZE_THREAD  "Enable ThreadSanitizer"                            OFF)
    option(LACE_SANITIZE_UB      "Enable UndefinedBehaviorSanitizer"                 OFF)
endif()

if(BUILD_SHARED_LIBS)
    message(WARNING
        "Building Lace as a shared library is not recommended. PLT indirection "
        "on hot paths hurts performance, and Lace must be instantiated exactly "
        "once per process. Instead, build Lace as a static library with LACE_ENABLE_PIC "
        "and link it into your shared library or executable."
    )
endif()

# ── Platform feature detection ───────────────────────────────────────────────

include(CheckSymbolExists)
check_symbol_exists(sched_getaffinity "sched.h" HAVE_SCHED_GETAFFINITY)
check_symbol_exists(mmap "sys/mman.h" HAVE_MMAP)

if(NOT WIN32 AND NOT HAVE_MMAP)
    message(FATAL_ERROR "mmap not found (required on non-Windows platforms)")
endif()

# ── Dependencies ─────────────────────────────────────────────────────────────

if(NOT MSVC)
    find_package(Threads REQUIRED)
endif()

if(LACE_USE_HWLOC)
    find_package(HWLOC REQUIRED)
endif()

# ── Sanitizer support (top-level builds only) ────────────────────────────────

if(lace_IS_TOP_LEVEL)
    if(LACE_SANITIZE_ADDRESS AND LACE_SANITIZE_THREAD)
        message(FATAL_ERROR "AddressSanitizer and ThreadSanitizer cannot be combined")
    endif()

    set(_san "")
    if(LACE_SANITIZE_ADDRESS)
        list(APPEND _san "address")
    endif()
    if(LACE_SANITIZE_THREAD)
        list(APPEND _san "thread")
    endif()
    if(LACE_SANITIZE_UB)
        list(APPEND _san "undefined")
    endif()

    if(_san)
        if(NOT CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
            message(FATAL_ERROR "Sanitizers are not supported with compiler '${CMAKE_C_COMPILER_ID}'")
        endif()
        string(REPLACE ";" "," _san_list "${_san}")
        set(_san_flags "-fsanitize=${_san_list} -fno-omit-frame-pointer")
        string(APPEND CMAKE_C_FLAGS             " ${_san_flags}")
        string(APPEND CMAKE_CXX_FLAGS           " ${_san_flags}")
        string(APPEND CMAKE_EXE_LINKER_FLAGS    " ${_san_flags}")
        string(APPEND CMAKE_SHARED_LINKER_FLAGS " ${_san_flags}")
        string(APPEND CMAKE_MODULE_LINKER_FLAGS " ${_san_flags}")
        message(STATUS "Sanitizers enabled: ${_san_list}")
    endif()
    unset(_san)
endif()

# ── Generated header (before targets, so FILE_SET can reference it) ──────────

configure_file(src/lace_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/lace_config.h @ONLY)

# ── Library targets ──────────────────────────────────────────────────────────

foreach(_target lace lace14)

    # ── Sources and public headers ───────────────────────────────────────
    add_library(${_target})
    add_library(lace::${_target} ALIAS ${_target})

    # The static header lives under src/ and the generated lace_config.h lives
    # under the build directory. These are kept in two separate file sets
    # (rather than one FILE_SET with two BASE_DIRS) because CMake requires
    # that all base directories accumulated for a single file set be disjoint
    # (no base directory may be a subdirectory of another); depending on how
    # the source and binary directories are laid out relative to each other,
    # combining them into one file set can trip that check.
    if(_target STREQUAL "lace")
        target_sources(${_target}
          PRIVATE src/lace.c
          PUBLIC  FILE_SET HEADERS
                  BASE_DIRS src
                  FILES src/lace.h
          PUBLIC  FILE_SET generated_headers
                  TYPE HEADERS
                  BASE_DIRS ${CMAKE_CURRENT_BINARY_DIR}
                  FILES ${CMAKE_CURRENT_BINARY_DIR}/lace_config.h
        )
    else()
        target_sources(${_target}
          PRIVATE src/lace14.c
          PUBLIC  FILE_SET HEADERS
                  BASE_DIRS src
                  FILES src/lace14.h
          PUBLIC  FILE_SET generated_headers
                  TYPE HEADERS
                  BASE_DIRS ${CMAKE_CURRENT_BINARY_DIR}
                  FILES ${CMAKE_CURRENT_BINARY_DIR}/lace_config.h
        )
    endif()

    # ── Compile features ─────────────────────────────────────────────────
    target_compile_features(${_target} PUBLIC c_std_11)

    # ── PIC ──────────────────────────────────────────────────────────────
    if(BUILD_SHARED_LIBS OR LACE_ENABLE_PIC)
        set_target_properties(${_target} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    endif()

    # ── Compiler flags ───────────────────────────────────────────────────

    # we're adding quite a few extra warnings, to ensure Lace as a library is compliant with the
    # C standard and won't trigger warnings or errors when included in a stricter project

    if(CMAKE_C_COMPILER_ID MATCHES "GNU")
        target_compile_options(${_target} PRIVATE
            -pipe
            -Wall -Wextra -Wpedantic -Werror
            -Wshadow # avoid potential shadowing bugs
            -Wnull-dereference # warn when we could null dereference
            -Wformat=2 # improve format string checking for logging
            -Wstrict-prototypes # functions without parameters must be void
            -Wmissing-prototypes # enforce non-static functions have prior declarations
            -Wconversion # catch implicit narrowing conversions
            -Wsign-conversion # catch implicit changes to signed/unsigned
            $<$<CONFIG:Debug>:-O0>
            $<$<BOOL:${LACE_NATIVE_OPT}>:-march=native>
        )
    elseif(CMAKE_C_COMPILER_ID MATCHES "Clang")
        target_compile_options(${_target} PRIVATE
            -pipe
            -Wall -Wextra -Wpedantic -Werror
            -Wshadow # avoid potential shadowing bugs
            -Wnull-dereference # warn when we could null dereference
            -Wformat=2 # improve format string checking for logging
            -Wstrict-prototypes # functions without parameters must be void
            -Wmissing-prototypes # enforce non-static functions have prior declarations
            -Wconversion # catch implicit narrowing conversions
            -Wsign-conversion # catch implicit changes to signed/unsigned
            -Watomic-implicit-seq-cst # ensure all atomic load/store are explicit
            $<$<CONFIG:Debug>:-O0>
            $<$<BOOL:${LACE_NATIVE_OPT}>:-march=native>
        )
    elseif(MSVC)
        target_compile_options(${_target} PUBLIC
            /Zc:__STDC__
            /Zc:preprocessor
            /experimental:c11atomics
        )
        target_compile_options(${_target} PRIVATE
            /W4
            /WX
            /wd4200 # Stop MSVC from complaining about standard C99/C11 code
        )
    endif()

    # ── Link dependencies ────────────────────────────────────────────────
    if(NOT MSVC)
        target_link_libraries(${_target} PUBLIC Threads::Threads)
    endif()

    if(WIN32)
        target_link_libraries(${_target} PUBLIC synchronization)
    endif()

    if(LACE_USE_HWLOC)
        target_link_libraries(${_target} PUBLIC HWLOC::hwloc)
    endif()

    # ── Shared library versioning ────────────────────────────────────────
    if(BUILD_SHARED_LIBS)
        set_target_properties(${_target} PROPERTIES
            VERSION   ${PROJECT_VERSION}
            SOVERSION ${PROJECT_VERSION_MAJOR}
            WINDOWS_EXPORT_ALL_SYMBOLS ON
        )
    endif()

endforeach()

# ── Tests and benchmarks (top-level only) ────────────────────────────────────

if(lace_IS_TOP_LEVEL)
    enable_language(CXX)

    if(LACE_BUILD_BENCHMARKS)
        add_subdirectory(benchmarks)
    endif()

    if(LACE_BUILD_TESTS)
        include(CTest)
        enable_testing()
        add_subdirectory(test)
    endif()
endif()

# ── Hide options when consumed as a subproject ───────────────────────────────

if(NOT lace_IS_TOP_LEVEL)
    mark_as_advanced(
        BUILD_SHARED_LIBS
        LACE_ENABLE_PIC
        LACE_NATIVE_OPT
        LACE_BACKOFF
        LACE_USE_HWLOC
        LACE_PIE_TIMES
        LACE_COUNT_TASKS
        LACE_COUNT_STEALS
        LACE_COUNT_SPLITS
    )
endif()

# ── Installation ─────────────────────────────────────────────────────────────

include(GNUInstallDirs)

install(TARGETS lace lace14
    EXPORT lace-targets
    ARCHIVE       DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY       DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME       DESTINATION ${CMAKE_INSTALL_BINDIR}
    FILE_SET HEADERS           DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILE_SET generated_headers DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT lace-targets
    FILE      lace-targets.cmake
    NAMESPACE lace::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lace
)

# ── CMake package config ─────────────────────────────────────────────────────

include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/lace-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/lace-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lace
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/lace-config-version.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMinorVersion
)

install(
    FILES
        ${CMAKE_CURRENT_BINARY_DIR}/cmake/lace-config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/cmake/lace-config-version.cmake
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindHWLOC.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lace
)

# ── pkg-config ───────────────────────────────────────────────────────────────

set(_PC_EXTRA_CFLAGS "")
set(_PC_EXTRA_LIBS "")

if(LACE_USE_HWLOC)
    foreach(inc ${HWLOC_INCLUDE_DIRS})
        list(APPEND _PC_EXTRA_CFLAGS "-I${inc}")
    endforeach()
    foreach(lib ${HWLOC_LIBRARIES})
        get_filename_component(libdir ${lib} DIRECTORY)
        get_filename_component(libname ${lib} NAME_WE)
        string(REGEX REPLACE "^lib" "" libname "${libname}")
        list(APPEND _PC_EXTRA_LIBS "-L${libdir}" "-l${libname}")
    endforeach()
endif()

if(MINGW)
    list(APPEND _PC_EXTRA_LIBS "-lsynchronization")
endif()

string(JOIN " " PC_EXTRA_CFLAGS ${_PC_EXTRA_CFLAGS})
string(JOIN " " PC_EXTRA_LIBS ${_PC_EXTRA_LIBS})

if(NOT MSVC)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/lace.pc.in"
                   "${CMAKE_CURRENT_BINARY_DIR}/lace.pc" @ONLY)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/lace14.pc.in"
                   "${CMAKE_CURRENT_BINARY_DIR}/lace14.pc" @ONLY)
    install(
        FILES
            ${CMAKE_CURRENT_BINARY_DIR}/lace.pc
            ${CMAKE_CURRENT_BINARY_DIR}/lace14.pc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
    )
endif()

# ── Summary ──────────────────────────────────────────────────────────────────

# # Summary
# message(STATUS "")
# message(STATUS "Lace ${PROJECT_VERSION} configuration summary:")
# message(STATUS "  Build type        : ${CMAKE_BUILD_TYPE}")
# message(STATUS "  Shared libraries  : ${BUILD_SHARED_LIBS}")
# message(STATUS "  Top-level project : ${lace_IS_TOP_LEVEL}")
# message(STATUS "  Use HWLOC         : ${LACE_USE_HWLOC}")
# message(STATUS "  Use -march=native : ${LACE_NATIVE_OPT}")
# message(STATUS "  Use -fPIC         : ${LACE_ENABLE_PIC}")
# message(STATUS "  Backoff when idle : ${LACE_BACKOFF}")
# message(STATUS "  Count tasks       : ${LACE_COUNT_TASKS}")
# message(STATUS "  Count steals      : ${LACE_COUNT_STEALS}")
# message(STATUS "  Count splits      : ${LACE_COUNT_SPLITS}")
# message(STATUS "  Record PIE times  : ${LACE_PIE_TIMES}")
# if(lace_is_top_level)
#     message(STATUS "  Build tests       : ${LACE_BUILD_TESTS}")
#     message(STATUS "  Build benchmarks  : ${LACE_BUILD_BENCHMARKS}")
#     message(STATUS "  Sanitize address  : ${LACE_SANITIZE_ADDRESS}")
#     message(STATUS "  Sanitize thread   : ${LACE_SANITIZE_THREAD}")
#     message(STATUS "  Sanitize UB       : ${LACE_SANITIZE_UB}")
# endif()
# message(STATUS "")
