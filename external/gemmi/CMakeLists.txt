cmake_minimum_required(VERSION 3.15...3.22)

# get version string from version.hpp
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/include/gemmi/version.hpp"
     gemmi_version_line REGEX "#define GEMMI_VERSION \"")
string(REGEX REPLACE ".+\"([0-9]+\.[0-9]+\.[0-9]+)(-dev)?\"" "\\1"
       gemmi_version_str ${gemmi_version_line})

project(gemmi LANGUAGES C CXX VERSION ${gemmi_version_str})
message(STATUS "Gemmi version ${PROJECT_VERSION}")

option(BUILD_SHARED_LIBS "Build shared library" ON) # currently, always ON
option(USE_ZLIB_NG "Use zlib-ng instead of zlib" OFF)
option(FETCH_ZLIB_NG "Download, build (statically) and use zlib-ng" OFF)
option(INTERNAL_ZLIB "Use subset of zlib distributed with gemmi" ON)

include(GNUInstallDirs)  # for CMAKE_INSTALL_LIBDIR, etc
include(CheckIncludeFile)  # for check_include_file

if (NOT DEFINED CMAKE_VISIBILITY_INLINES_HIDDEN)
  set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)
endif()
if (NOT DEFINED CMAKE_CXX_VISIBILITY_PRESET)
  set(CMAKE_CXX_VISIBILITY_PRESET hidden)
endif()

# Unless it's explicitly specifed, we want LTO for non-debug builds only.
if (NOT DEFINED CMAKE_INTERPROCEDURAL_OPTIMIZATION)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT is_ipo_supported OUTPUT ipo_error)
  if(is_ipo_supported)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_MINSIZEREL ON)
  else()
    message(STATUS "Interprocedural optimization not supported: ${error}")
  endif()
endif()

if (DEFINED ENV{CXXFLAGS})
  set(USING_ENV_CXXFLAGS ON CACHE BOOL "" FORCE)
endif()

# always build in Release mode
set(CMAKE_BUILD_TYPE Release CACHE STRING "Force Release build" FORCE)

if (FETCH_ZLIB_NG)
  message(STATUS " =============== fetching and building zlib-ng ===============")
  include(FetchContent)
  FetchContent_Declare(zlib-ng GIT_REPOSITORY https://github.com/zlib-ng/zlib-ng.git
                               GIT_TAG        2.1.x
                               GIT_SHALLOW TRUE)
  FetchContent_GetProperties(zlib-ng)
  if (NOT zlib-ng_POPULATED)
    FetchContent_Populate(zlib-ng)
    set(ZLIB_ENABLE_TESTS OFF CACHE INTERNAL "")
    set(BUILD_SHARED_LIBS_here ${BUILD_SHARED_LIBS})
    set(BUILD_SHARED_LIBS OFF CACHE INTERNAL "build static zlib-ng")
    add_subdirectory(${zlib-ng_SOURCE_DIR} ${zlib-ng_BINARY_DIR} EXCLUDE_FROM_ALL)
    set(BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS_here} CACHE BOOL "Build shared library" FORCE)
    set_property(TARGET zlib PROPERTY POSITION_INDEPENDENT_CODE ON)
    # zlib-ng, unless built WITH_NATIVE_INSTRUCTIONS=ON (-march=native),
    # builds multiple variants of some functions, using different CPU features.
    # IPO (LTO) has been reported to mess it up leading to illegal-instruction
    # crash. When zlib-ng is built as a standalone project, not subdirectory,
    # these flags are set in zlib-ng/CMakeLists.txt.
    if (NOT WITH_NATIVE_INSTRUCTIONS)
      set_target_properties(zlib PROPERTIES INTERPROCEDURAL_OPTIMIZATION OFF
                                            INTERPROCEDURAL_OPTIMIZATION_RELEASE OFF
                                            INTERPROCEDURAL_OPTIMIZATION_MINSIZEREL OFF)
    endif()
    message(STATUS " =============== end of zlib-ng configuration ===============\n")
  endif()
elseif (USE_ZLIB_NG)
  set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/tools")
  find_package(ZLIB_NG REQUIRED)
else()
  if (NOT INTERNAL_ZLIB)
    # CMake >= 3.24 respects variable ZLIB_USE_STATIC_LIBS
    find_package(ZLIB)
  endif()
  if (NOT ZLIB_FOUND)
    message(STATUS "The build will use zlib code from third_party/zlib.")
  endif()
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

add_library(gemmi_headers INTERFACE)
add_library(gemmi::headers ALIAS gemmi_headers)
target_include_directories(gemmi_headers INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
target_compile_features(gemmi_headers INTERFACE cxx_std_11)
set_target_properties(gemmi_headers PROPERTIES EXPORT_NAME headers)

add_library(gemmi_cpp SHARED
            src/ddl.cpp src/gz.cpp src/json.cpp src/mmcif.cpp src/mmread_gz.cpp
            src/polyheur.cpp src/read_cif.cpp src/resinfo.cpp)

add_library(gemmi::gemmi_cpp ALIAS gemmi_cpp)
set_property(TARGET gemmi_cpp PROPERTY POSITION_INDEPENDENT_CODE ON)
# set_property(TARGET gemmi_cpp PROPERTY CXX_VISIBILITY_PRESET hidden)
# target_compile_definitions(gemmi_cpp PUBLIC GEMMI_BUILD)
target_include_directories(gemmi_cpp PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/zlib
)

target_compile_definitions(gemmi_cpp PUBLIC GEMMI_SHARED)
target_link_libraries(gemmi_cpp PUBLIC gemmi_headers)

if (FETCH_ZLIB_NG)
  target_compile_definitions(gemmi_cpp PUBLIC USE_ZLIB_NG=1)
  target_link_libraries(gemmi_cpp PUBLIC zlib)
elseif (ZLIB_NG_FOUND)
  target_compile_definitions(gemmi_cpp PUBLIC USE_ZLIB_NG=1)
  target_link_libraries(gemmi_cpp PUBLIC ZLIB_NG::ZLIB_NG)
elseif (ZLIB_FOUND)
  target_link_libraries(gemmi_cpp PUBLIC ZLIB::ZLIB)
else()
  message(STATUS "->: Using internal zlib code.")
  target_sources(gemmi_cpp PUBLIC
      "third_party/zlib/adler32.c"
      "third_party/zlib/crc32.c"
      "third_party/zlib/gzlib.c"
      "third_party/zlib/gzread.c"
      "third_party/zlib/inflate.c"
      "third_party/zlib/inftrees.c"
      "third_party/zlib/inffast.c"
      "third_party/zlib/zutil.c")
  check_include_file(unistd.h has_unistd_h)
  target_compile_definitions(gemmi_cpp PUBLIC NO_GZCOMPRESS=1 DYNAMIC_CRC_TABLE=1)
  if (has_unistd_h)
    target_compile_definitions(gemmi_cpp PUBLIC Z_HAVE_UNISTD_H=1)
  endif()
  target_include_directories(gemmi_cpp PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/third_party/zlib")
  if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    # /wd4267 - conversion from 'size_t' to 'unsigned int', possible loss of data
    target_compile_options(gemmi_cpp PUBLIC /wd4267)
  endif()
endif()

