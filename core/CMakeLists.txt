cmake_minimum_required(VERSION 3.15...3.22)
project(lahuta VERSION 1.0.0 LANGUAGES CXX C)

include(GNUInstallDirs)
include(CheckIncludeFile)

set(HWY_ENABLE_TESTS      OFF CACHE BOOL "Enable HWY tests" FORCE)
set(HWY_ENABLE_EXAMPLES   OFF CACHE BOOL "Disable highway examples" FORCE)
set(HWY_ENABLE_BENCHMARKS OFF CACHE BOOL "Disable highway benchmarks" FORCE)
set(HWY_ENABLE_CONTRIB    OFF CACHE BOOL "Disable highway contrib targets" FORCE)
set(HWY_ENABLE_INSTALL    OFF CACHE BOOL "Disable highway install targets" FORCE)

set(SPDLOG_FMT_EXTERNAL      OFF CACHE BOOL "Use bundled fmt")
set(SPDLOG_BUILD_EXAMPLES    OFF CACHE BOOL "")
set(SPDLOG_BUILD_TESTS       OFF CACHE BOOL "")
set(SPDLOG_BUILD_EXTRAS      OFF CACHE BOOL "")
set(SPDLOG_BUILD_PERF_TESTS  OFF CACHE BOOL "")
set(SPDLOG_BUILD_HEADER_ONLY ON  CACHE BOOL "")

option(USE_ZLIB_NG    "Use zlib-ng instead of zlib" OFF)
option(FETCH_ZLIB_NG  "Download, build (statically) and use zlib-ng" OFF)
option(INTERNAL_ZLIB  "Use subset of zlib distributed with gemmi" OFF)
option(GENERATE_STUBS "Generate Python type stubs" OFF)
option(USE_EIGEN3     "Do not use Eigen3" ON)

set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_BUILD_RPATH_USE_ORIGIN ON)

if (DEFINED ENV{CXXFLAGS})
  set(USING_ENV_CXXFLAGS ON CACHE BOOL "" FORCE)
endif()

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release)
endif()

# We need to find zlib before adding gemmi so it uses the same zlib
find_package(ZLIB REQUIRED)
if(ZLIB_FOUND)
  message(STATUS "Found zlib: ${ZLIB_LIBRARIES}")
else()
  message(FATAL_ERROR "zlib not found - required for compression")
endif()

if(DEFINED ENV{BOOST_ROOT} AND NOT DEFINED BOOST_ROOT)
  set(BOOST_ROOT $ENV{BOOST_ROOT})
endif()

set(_BOOST_HINTS ${BOOST_ROOT} $ENV{CONDA_PREFIX} /usr/local /opt/homebrew /usr)

find_package(Boost 1.81 REQUIRED COMPONENTS graph HINTS ${_BOOST_HINTS} PATH_SUFFIXES lib64 lib)

if(NOT Boost_FOUND)
  message(FATAL_ERROR "Could not find Boost≥1.81."
  "Tried:\n -DBOOST_ROOT=${BOOST_ROOT}\n $ENV{CONDA_PREFIX}\n /usr/local, /opt/homebrew, /usr\n"
  "Set BOOST_ROOT to the correct path or install Boost.")
endif()

ADD_LIBRARY(lmdb STATIC
    ${CMAKE_SOURCE_DIR}/external/lmdb/mdb.c
    ${CMAKE_SOURCE_DIR}/external/lmdb/midl.c
)

set_target_properties(lmdb PROPERTIES LINKER_LANGUAGE C)
TARGET_INCLUDE_DIRECTORIES(lmdb PUBLIC ${CMAKE_SOURCE_DIR}/external)

add_subdirectory(external/spdlog)
# gemmi first to inherit its zlib configuration
add_subdirectory(external/gemmi)
add_subdirectory(external/rdkit)
add_subdirectory(external/foldseek)
add_dependencies(foldseek gemmi_cpp)
add_subdirectory(external/distopia)

# set(CMAKE_BUILD_TYPE Debug)
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

add_library(lahuta_cpp STATIC
            src/bonds.cpp
            src/bond_order.cpp
            src/convert.cpp
            src/nsgrid.cpp
            src/lahuta.cpp
            src/visitor.cpp
            src/struct_unit.cpp
            src/topology.cpp
            src/residues.cpp
            src/aromatics.cpp
            src/luni_props.cpp
            src/logging.cpp

            src/models/topology.cpp
            src/models/parser.cpp
            src/models/topology_impl.cpp
            src/models/topology/engine.cpp
            src/models/topology/kernels/atoms.cpp
            src/models/topology/kernels/bonds.cpp
            src/models/topology/kernels/positions.cpp
            src/models/topology/kernels/aromatics.cpp
            src/models/topology/kernels/disulfides.cpp
            src/models/topology/kernels/build.cpp

            external/mmap/MemoryMapped.cpp
            src/valence_model.cpp

            src/ob/kekulize.cpp
            src/ob/clean_mol.cpp
            src/ob/bitvec.cpp

            src/contacts/aromaticity.cpp

            src/chemistry/group_typing.cpp
            src/chemistry/types/charges.cpp
            src/chemistry/types/metalics.cpp
            src/chemistry/types/hydrophobic.cpp
            src/chemistry/types/halogens.cpp
            src/chemistry/types/hbonding.cpp

            src/mapping/bitset_mapping.cpp
            src/mapping/mapping_matrix.cpp
            src/mapping/residue_set.cpp

            src/topology/engine.cpp
            src/topology/kernels/neighbors.cpp
            src/topology/kernels/bonds.cpp
            src/topology/kernels/rings.cpp
            src/topology/kernels/residues.cpp
            src/topology/kernels/atom_typing.cpp

            src/chemistry/neighbors.cpp
            src/chemistry/geometry.cpp
            src/chemistry/predicates.cpp

            src/entities/contact.cpp

            src/contacts/molstar/ionic.cpp
            src/contacts/molstar/hbond.cpp
            src/contacts/molstar/hydrophobic.cpp
            src/contacts/molstar/halogen.cpp
            src/contacts/molstar/cationpi.cpp
            src/contacts/molstar/pistacking.cpp
            src/contacts/molstar/metalic.cpp

            src/contacts/arpeggio/hbonds.cpp
            src/contacts/arpeggio/hydrophobic.cpp
            src/contacts/arpeggio/ionic.cpp
            src/contacts/arpeggio/aromatic.cpp
            src/contacts/arpeggio/carbonyl.cpp
            src/contacts/arpeggio/atom_plane.cpp

            src/cli/arg_validation.cpp
            src/cli/global_options.cpp
            src/commands/contacts.cpp
            src/commands/createdb.cpp
)

set_property(TARGET lahuta_cpp PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(lahuta_cpp PUBLIC gemmi_cpp foldseek rdkit libdistopia lmdb ZLIB::ZLIB spdlog::spdlog_header_only Boost::headers)

target_include_directories(lahuta_cpp PUBLIC
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/external
  ${CMAKE_SOURCE_DIR}/external/foldseek
  ${CMAKE_SOURCE_DIR}/external/foldseek/generated
  ${CMAKE_SOURCE_DIR}/external/gemmi/include
  ${CMAKE_SOURCE_DIR}/external/gemmi/third_party
  ${CMAKE_SOURCE_DIR}/external/gemmi/third_party/zlib
  ${CMAKE_SOURCE_DIR}/external/distopia/include
  ${CMAKE_SOURCE_DIR}/external/lmdb
)

add_library(lahuta_shared_headers INTERFACE)
target_include_directories(lahuta_shared_headers INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/src/fseek
  ${CMAKE_CURRENT_SOURCE_DIR}/external/ctpl
  ${CMAKE_CURRENT_SOURCE_DIR}/examples
)

find_package(Eigen3 QUIET)
if(USE_EIGEN3)
  if(NOT EIGEN3_FOUND)
    add_subdirectory(external/Eigen)
    add_definitions(-DLAHUTA_USE_EIGEN=1)
  endif()
endif()

if(TARGET Eigen3::Eigen)
  message(STATUS "Eigen3 found.")
  target_link_libraries(lahuta_cpp PUBLIC Eigen3::Eigen)
  set(LAHUTA_HAS_EIGEN ON)
endif()

function(add_lahuta_executable target main_file)
  add_executable(${target} ${main_file})

  if (LAHUTA_HAS_EIGEN)
    target_link_libraries(${target} PRIVATE Eigen3::Eigen)
  endif()

  target_link_libraries(${target} PUBLIC lahuta_cpp ZLIB::ZLIB spdlog::spdlog_header_only)
endfunction()

add_lahuta_executable(lahuta apps/main.cpp)
