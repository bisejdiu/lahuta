cmake_minimum_required(VERSION 3.20)
project(lahuta VERSION 1.0.0 LANGUAGES CXX C)

include(CTest) # Provides standard BUILD_TESTING option
include(GNUInstallDirs)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(LAHUTA_BUILD_CLI            "Build CLI executables" ON)
option(LAHUTA_BUILD_PYTHON         "Build Python bindings" ON)
option(LAHUTA_GENERATE_PY_STUBS    "Generate Python type stubs with pybind11_stubgen" OFF)
option(LAHUTA_WITH_EIGEN           "Use Eigen3 (system or bundled fallback)" ON)

set(SPDLOG_FMT_EXTERNAL      OFF CACHE BOOL "Use bundled fmt for spdlog")
set(SPDLOG_BUILD_EXAMPLES    OFF CACHE BOOL "" )
set(SPDLOG_BUILD_TESTS       OFF CACHE BOOL "" )
set(SPDLOG_BUILD_EXTRAS      OFF CACHE BOOL "" )
set(SPDLOG_BUILD_PERF_TESTS  OFF CACHE BOOL "" )
set(SPDLOG_BUILD_HEADER_ONLY ON  CACHE BOOL "" )

# Highway (used by distopia). Ensure it does not build extras/tests.
set(HWY_ENABLE_TESTS      OFF CACHE BOOL "Disable highway tests"           FORCE)
set(HWY_ENABLE_EXAMPLES   OFF CACHE BOOL "Disable highway examples"        FORCE)
set(HWY_ENABLE_BENCHMARKS OFF CACHE BOOL "Disable highway benchmarks"      FORCE)
set(HWY_ENABLE_CONTRIB    OFF CACHE BOOL "Disable highway contrib targets" FORCE)
set(HWY_ENABLE_INSTALL    OFF CACHE BOOL "Disable highway install targets" FORCE)

find_package(ZLIB REQUIRED)

# Boost search hints
if(DEFINED ENV{BOOST_ROOT} AND NOT DEFINED BOOST_ROOT)
  set(BOOST_ROOT $ENV{BOOST_ROOT})
endif()
set(_BOOST_HINTS ${BOOST_ROOT} $ENV{CONDA_PREFIX} /usr/local /opt/homebrew /usr)
find_package(Boost 1.81 REQUIRED COMPONENTS graph HINTS ${_BOOST_HINTS} PATH_SUFFIXES lib64 lib)

# Interface header collection used by foldseek
add_library(lahuta_shared_headers INTERFACE)
target_include_directories(lahuta_shared_headers INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/src/fseek
  ${CMAKE_CURRENT_SOURCE_DIR}/external/ctpl
  ${CMAKE_CURRENT_SOURCE_DIR}/examples
)

add_subdirectory(external/spdlog)
add_subdirectory(external/gemmi)    # Build gemmi first (zlib alignment)
add_subdirectory(external/rdkit)
add_subdirectory(external/foldseek)
add_dependencies(foldseek gemmi_cpp)
add_subdirectory(external/distopia)

# LMDB
add_library(lmdb STATIC
  ${CMAKE_CURRENT_SOURCE_DIR}/external/lmdb/mdb.c
  ${CMAKE_CURRENT_SOURCE_DIR}/external/lmdb/midl.c
)
set_target_properties(lmdb PROPERTIES LINKER_LANGUAGE C)
target_include_directories(lmdb PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/external)

# Core library
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

target_compile_features(lahuta_core PUBLIC cxx_std_17)
set_property(TARGET lahuta_core PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(lahuta_core PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/external
  ${CMAKE_CURRENT_SOURCE_DIR}/external/foldseek
  ${CMAKE_CURRENT_SOURCE_DIR}/external/foldseek/generated
  ${CMAKE_CURRENT_SOURCE_DIR}/external/gemmi/include
  ${CMAKE_CURRENT_SOURCE_DIR}/external/gemmi/third_party
  ${CMAKE_CURRENT_SOURCE_DIR}/external/gemmi/third_party/zlib
  ${CMAKE_CURRENT_SOURCE_DIR}/external/distopia/include
  ${CMAKE_CURRENT_SOURCE_DIR}/external/lmdb
)

target_link_libraries(lahuta_core PUBLIC
  gemmi_cpp
  foldseek
  rdkit
  libdistopia
  lmdb
  ZLIB::ZLIB
  spdlog::spdlog_header_only
  Boost::headers
)

# Eigen (system first, fallback to bundled). Define LAHUTA_USE_EIGEN=1 when linked.
if(LAHUTA_WITH_EIGEN)
  find_package(Eigen3 QUIET)
  if(NOT TARGET Eigen3::Eigen)
    add_subdirectory(external/Eigen)
  endif()
  if(TARGET Eigen3::Eigen)
    target_link_libraries(lahuta_core PUBLIC Eigen3::Eigen)
    target_compile_definitions(lahuta_core PUBLIC LAHUTA_USE_EIGEN=1)
    set(LAHUTA_HAS_EIGEN ON)
  endif()
endif()

# CLI executables
function(add_lahuta_executable target main_file)
  add_executable(${target} ${main_file})
  target_compile_features(${target} PRIVATE cxx_std_17)
  target_link_libraries(${target} PRIVATE lahuta_core ZLIB::ZLIB spdlog::spdlog_header_only)
endfunction()

if(LAHUTA_BUILD_CLI)
  add_lahuta_executable(lahuta_cli apps/main.cpp)
  set_target_properties(lahuta_cli PROPERTIES OUTPUT_NAME "lahuta")
  add_lahuta_executable(resolver apps/use_resolver.cpp)
endif()

# Python bindings
if(LAHUTA_BUILD_PYTHON)
  find_package(pybind11 QUIET CONFIG)
  if(pybind11_FOUND)
    message(STATUS "pybind11 found. Building Python bindings.")
    find_package(Python3 COMPONENTS Interpreter Development.Module QUIET)

    set(PY_SRC_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/src)
    set(PY_PKG_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python)
    set(PY_PKG_LIB  ${PY_PKG_ROOT}/lahuta/lib)

    set(LAHUTA_PY_SOURCES
      ${PY_SRC_DIR}/lahuta.cpp

      ${PY_SRC_DIR}/utilities.cpp
      ${PY_SRC_DIR}/logger.cpp

      ${PY_SRC_DIR}/luni.cpp
      ${PY_SRC_DIR}/topology.cpp
      ${PY_SRC_DIR}/entity_id.cpp
      ${PY_SRC_DIR}/entities.cpp
      ${PY_SRC_DIR}/atom_types.cpp
      ${PY_SRC_DIR}/records.cpp
      ${PY_SRC_DIR}/properties.cpp

      ${PY_SRC_DIR}/grid_search.cpp
      ${PY_SRC_DIR}/distance.cpp
      ${PY_SRC_DIR}/contacts.cpp

      ${PY_SRC_DIR}/rdkit.cpp

      ${PY_SRC_DIR}/pipeline_dynamic.cpp
    )

    pybind11_add_module(lahuta_ext MODULE ${LAHUTA_PY_SOURCES})
    target_include_directories(lahuta_ext PRIVATE ${PY_SRC_DIR})
    target_link_libraries(lahuta_ext PRIVATE lahuta_core spdlog::spdlog_header_only ZLIB::ZLIB)
    target_compile_features(lahuta_ext PRIVATE cxx_std_17)
    set_target_properties(lahuta_ext PROPERTIES OUTPUT_NAME "lahuta")
    # Portable runtime search paths for loading dependent libs from the package tree
    set_target_properties(lahuta_ext PROPERTIES
      BUILD_RPATH "$<$<PLATFORM_ID:Darwin>:@loader_path;@loader_path/..>$<$<PLATFORM_ID:Linux>:$ORIGIN;$ORIGIN/..>"
    )

    # mapping submodule, depends at Python-import level on lahuta
    set(LAHUTA_MAPPING_SOURCES
      ${PY_SRC_DIR}/mapping/mapping.cpp
      ${PY_SRC_DIR}/mapping/fseek.cpp
      ${PY_SRC_DIR}/mapping/align.cpp
      ${PY_SRC_DIR}/mapping/te.cpp
    )
    pybind11_add_module(mapping_ext MODULE ${LAHUTA_MAPPING_SOURCES})
    target_include_directories(mapping_ext PRIVATE ${PY_SRC_DIR})
    target_link_libraries(mapping_ext PRIVATE lahuta_core spdlog::spdlog_header_only ZLIB::ZLIB)
    target_compile_features(mapping_ext PRIVATE cxx_std_17)
    set_target_properties(mapping_ext PROPERTIES OUTPUT_NAME "mapping")
    set_target_properties(mapping_ext PROPERTIES
      BUILD_RPATH "$<$<PLATFORM_ID:Darwin>:@loader_path;@loader_path/..>$<$<PLATFORM_ID:Linux>:$ORIGIN;$ORIGIN/..>"
    )

    # Copy extensions into the Python package tree, as individual custom targets
    add_custom_target(copy_lahuta_ext ALL
      COMMAND ${CMAKE_COMMAND} -E make_directory ${PY_PKG_LIB}
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:lahuta_ext> ${PY_PKG_LIB}/
      DEPENDS lahuta_ext
      COMMENT "Copying lahuta extension to package lib directory"
    )
    add_custom_target(copy_mapping_ext ALL
      COMMAND ${CMAKE_COMMAND} -E make_directory ${PY_PKG_LIB}
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:mapping_ext> ${PY_PKG_LIB}/
      DEPENDS mapping_ext copy_lahuta_ext
      COMMENT "Copying mapping extension to package lib directory"
    )

    # Stub generation (custom target)
    if(LAHUTA_GENERATE_PY_STUBS)
      if(NOT Python3_Interpreter_FOUND)
        message(WARNING "Python3 interpreter not found. Cannot generate type stubs for Python bindings.")
      else()
        add_custom_target(python_stubs ALL
          COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${PY_PKG_ROOT}
                  ${Python3_EXECUTABLE} -m pybind11_stubgen lahuta.lib.lahuta
                  --output-dir ${PY_PKG_ROOT}
                  --root-suffix=
                  --ignore-all-errors
                  --numpy-array-use-type-var
          COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${PY_PKG_ROOT}
                  ${Python3_EXECUTABLE} -m pybind11_stubgen lahuta.lib.mapping
                  --output-dir ${PY_PKG_ROOT}
                  --root-suffix=
                  --ignore-all-errors
                  --numpy-array-use-type-var
          WORKING_DIRECTORY ${PY_PKG_ROOT}
          DEPENDS copy_lahuta_ext copy_mapping_ext
          COMMENT "Generating type stubs for lahuta modules"
        )
      endif()
    endif()
  else()
    message(STATUS "pybind11 not found. Skipping Python bindings build.")
  endif()
endif()

# Tests
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
