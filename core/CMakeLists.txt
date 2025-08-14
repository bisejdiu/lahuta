cmake_minimum_required(VERSION 3.15...3.22)
project(lahuta VERSION 1.0.0 LANGUAGES CXX C)
include(CTest)

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
option(GENERATE_STUBS "Generate Python type stubs" ON)
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
  message(FATAL_ERROR "Could not find Boost>=1.81."
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

add_library(lahuta_core STATIC
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
            src/topology/kernels/seed_from_model.cpp

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

set_property(TARGET lahuta_core PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(lahuta_core PUBLIC gemmi_cpp foldseek rdkit libdistopia lmdb ZLIB::ZLIB spdlog::spdlog_header_only Boost::headers)

target_include_directories(lahuta_core PUBLIC
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
  target_link_libraries(lahuta_core PUBLIC Eigen3::Eigen)
  set(LAHUTA_HAS_EIGEN ON)
endif()

function(add_lahuta_executable target main_file)
  add_executable(${target} ${main_file})

  if (LAHUTA_HAS_EIGEN)
    target_link_libraries(${target} PRIVATE Eigen3::Eigen)
  endif()

  target_link_libraries(${target} PUBLIC lahuta_core ZLIB::ZLIB spdlog::spdlog_header_only)
endfunction()

# CLI executable
add_lahuta_executable(lahuta_cli apps/main.cpp)
add_lahuta_executable(resolver apps/use_resolver.cpp)
set_target_properties(lahuta_cli PROPERTIES OUTPUT_NAME "lahuta")

find_package(pybind11 QUIET CONFIG)
if(pybind11_FOUND)
  message(STATUS "pybind11 found. Building Python bindings.")

  find_package(Python3 COMPONENTS Interpreter Development.Module QUIET)

  set(LAHUTA_PY_SOURCES
      ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/src/lahuta.cpp

      ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/src/utilities.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/src/logger.cpp

      ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/src/luni.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/src/topology.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/src/entity_id.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/src/entities.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/src/atom_types.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/src/records.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/src/properties.cpp

      ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/src/grid_search.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/src/distance.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/src/contacts.cpp

      ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/src/rdkit.cpp

      ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/src/pipeline_dynamic.cpp
  )

  # Python module
  pybind11_add_module(lahuta_ext MODULE ${LAHUTA_PY_SOURCES})

  target_include_directories(lahuta_ext PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/src)
  target_link_libraries(lahuta_ext PRIVATE lahuta_core spdlog::spdlog_header_only ZLIB::ZLIB)

  set_target_properties(lahuta_ext PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
  set_target_properties(lahuta_ext PROPERTIES BUILD_RPATH "@loader_path;@loader_path/..")
  set_target_properties(lahuta_ext PROPERTIES OUTPUT_NAME "lahuta")

  add_custom_command(TARGET lahuta_ext POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:lahuta_ext> 
            ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/lahuta/lib/
    COMMENT "Copying Python extension to package lib directory"
  )

  # Build mapping submodule
  set(LAHUTA_MAPPING_SOURCES
      ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/src/mapping/mapping.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/src/mapping/fseek.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/src/mapping/align.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/src/mapping/te.cpp
  )

  pybind11_add_module(mapping_ext MODULE ${LAHUTA_MAPPING_SOURCES})

  target_include_directories(mapping_ext PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/src)
  target_link_libraries(mapping_ext PRIVATE lahuta_core spdlog::spdlog_header_only ZLIB::ZLIB)

  set_target_properties(mapping_ext PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
  set_target_properties(mapping_ext PROPERTIES BUILD_RPATH "@loader_path;@loader_path/..")
  set_target_properties(mapping_ext PROPERTIES OUTPUT_NAME "mapping")

  add_custom_command(TARGET mapping_ext POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:mapping_ext>
            ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/lahuta/lib/
    COMMENT "Copying mapping extension to package lib directory"
  )

  if(GENERATE_STUBS)
    if(NOT Python3_Interpreter_FOUND)
      message(WARNING "Python3 interpreter not found. Cannot generate type stubs for Python bindings.")
    else()
      # Generate .pyi stubs for the built module next to the copied extension
      # We import the module as 'lahuta.lib.lahuta' by pointing PYTHONPATH to the interop/python root
      add_custom_command(TARGET lahuta_ext POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E env
                PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../interop/python
                ${Python3_EXECUTABLE} -m pybind11_stubgen lahuta.lib.lahuta
                --output-dir ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python
                --root-suffix=
                --ignore-all-errors
                --numpy-array-use-type-var
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python
        COMMENT "Generating type stubs for lahuta.lib.lahuta"
      )

      # Generate .pyi stubs for the mapping module
      add_custom_command(TARGET mapping_ext POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E env
                PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/../interop/python
                ${Python3_EXECUTABLE} -m pybind11_stubgen lahuta.lib.mapping
                --output-dir ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python
                --root-suffix=
                --ignore-all-errors
                --numpy-array-use-type-var
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python
        COMMENT "Generating type stubs for lahuta.lib.mapping"
      )
    endif()
  endif()

else()
  message(STATUS "pybind11 not found. Skipping Python bindings build.")
endif()

# Tests
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
