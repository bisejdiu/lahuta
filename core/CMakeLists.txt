cmake_minimum_required(VERSION 3.20)
project(lahuta VERSION 1.0.0 LANGUAGES CXX C)

include(CTest) # Provides BUILD_TESTING option
include(GNUInstallDirs)
include(CheckCXXSourceCompiles)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_CONFIGURATION_TYPES)
  set(_valid_build_types Debug Release RelWithDebInfo MinSizeRel)
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  endif()
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${_valid_build_types})
endif()

set(LAHUTA_OPT_LEVEL 3 CACHE STRING "Optimization level for Release/RelWithDebInfo (0/1/2/3 or 'fast')")
set_property(CACHE LAHUTA_OPT_LEVEL PROPERTY STRINGS 0 1 2 3 fast)
option(LAHUTA_ENABLE_LTO "Enable interprocedural optimization (LTO/IPO) for non-Debug builds" ON)
option(LAHUTA_NATIVE_ARCH "Enable -march=native for non-Debug builds" OFF)

# Apply interprocedural optimization consistently to subprojects
if(LAHUTA_ENABLE_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT _ipo_supported OUTPUT _ipo_error)
  if(_ipo_supported)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO ON)
  else()
    message(STATUS "IPO/LTO not supported by toolchain: ${_ipo_error}")
  endif()
endif()

if(NOT MSVC)
  set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g3 -fno-omit-frame-pointer")
  set(CMAKE_C_FLAGS_DEBUG     "${CMAKE_C_FLAGS_DEBUG}   -O0 -g3 -fno-omit-frame-pointer")

  set(_LAHUTA_OFLAG "")
  if(LAHUTA_OPT_LEVEL STREQUAL "fast")
    set(_LAHUTA_OFLAG "-Ofast")
  else()
    set(_LAHUTA_OFLAG "-O${LAHUTA_OPT_LEVEL}")
  endif()

  set(CMAKE_CXX_FLAGS_RELEASE         "${CMAKE_CXX_FLAGS_RELEASE}  ${_LAHUTA_OFLAG} -DNDEBUG -fno-omit-frame-pointer")
  set(CMAKE_C_FLAGS_RELEASE           "${CMAKE_C_FLAGS_RELEASE}    ${_LAHUTA_OFLAG} -DNDEBUG -fno-omit-frame-pointer")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${_LAHUTA_OFLAG} -g -fno-omit-frame-pointer")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO    "${CMAKE_C_FLAGS_RELWITHDEBINFO}   ${_LAHUTA_OFLAG} -g -fno-omit-frame-pointer")

  if(LAHUTA_NATIVE_ARCH)
    set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native")
    set(CMAKE_C_FLAGS_RELEASE          "${CMAKE_C_FLAGS_RELEASE}   -march=native -mtune=native")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -march=native -mtune=native")
    set(CMAKE_C_FLAGS_RELWITHDEBINFO   "${CMAKE_C_FLAGS_RELWITHDEBINFO}   -march=native -mtune=native")
  endif()
endif()

option(LAHUTA_BUILD_CLI            "Build CLI executables" ON)
option(LAHUTA_BUILD_PYTHON         "Build Python bindings" ON)
option(LAHUTA_GENERATE_PY_STUBS    "Generate Python type stubs with pybind11_stubgen" OFF)
option(LAHUTA_WITH_EIGEN           "Use Eigen3 (system or bundled fallback)" OFF)
option(LAHUTA_VENDOR_BOOST         "Use vendored Boost headers instead of system Boost" ON)

set(RDK_USE_BOOST_SERIALIZATION OFF CACHE BOOL "Disable rdkit boost serialization" FORCE)

set(SPDLOG_FMT_EXTERNAL      OFF CACHE BOOL "Use bundled fmt for spdlog")
set(SPDLOG_BUILD_EXAMPLES    OFF CACHE BOOL "" )
set(SPDLOG_BUILD_TESTS       OFF CACHE BOOL "" )
set(SPDLOG_BUILD_EXTRAS      OFF CACHE BOOL "" )
set(SPDLOG_BUILD_PERF_TESTS  OFF CACHE BOOL "" )
set(SPDLOG_BUILD_HEADER_ONLY ON  CACHE BOOL "" )

# Highway (used by distopia). Do not build extras/tests.
set(HWY_ENABLE_TESTS      OFF CACHE BOOL "Disable highway tests"           FORCE)
set(HWY_ENABLE_EXAMPLES   OFF CACHE BOOL "Disable highway examples"        FORCE)
set(HWY_ENABLE_BENCHMARKS OFF CACHE BOOL "Disable highway benchmarks"      FORCE)
set(HWY_ENABLE_CONTRIB    OFF CACHE BOOL "Disable highway contrib targets" FORCE)
set(HWY_ENABLE_INSTALL    OFF CACHE BOOL "Disable highway install targets" FORCE)

option(LAHUTA_DISTOPIA_ENABLE_SVE "Enable Highway SVE targets for distopia" OFF)

find_package(ZLIB REQUIRED)

# Boost headers
if(LAHUTA_VENDOR_BOOST)
  set(LAHUTA_BOOST_VENDOR_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/boost")
  if(NOT EXISTS "${LAHUTA_BOOST_VENDOR_DIR}/boost/version.hpp")
    message(FATAL_ERROR "Vendored Boost headers not found at ${CMAKE_CURRENT_SOURCE_DIR}/external/boost.\n"
                        "Ensure submodules are initialized (git submodule update --init --recursive),\n"
                        "or disable vendoring with -DLAHUTA_VENDOR_BOOST=OFF to use system Boost.")
  endif()
  add_library(Boost::headers INTERFACE IMPORTED)
  set_target_properties(Boost::headers PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${LAHUTA_BOOST_VENDOR_DIR}")
  message(STATUS "Using vendored Boost headers at ${LAHUTA_BOOST_VENDOR_DIR}")
else()
  if(DEFINED ENV{BOOST_ROOT} AND NOT DEFINED BOOST_ROOT)
    set(BOOST_ROOT $ENV{BOOST_ROOT})
  endif()
  set(_BOOST_HINTS ${BOOST_ROOT} $ENV{CONDA_PREFIX} /usr/local /opt/homebrew /usr)
  find_package(Boost 1.81 REQUIRED COMPONENTS graph HINTS ${_BOOST_HINTS} PATH_SUFFIXES lib64 lib)
endif()

#
# Quick check for Boost Graph headers, as without them lahuta cannot build.
#
set(_LAHUTA_BGL_INCLUDE_DIRS "")
if(LAHUTA_VENDOR_BOOST)
  list(APPEND _LAHUTA_BGL_INCLUDE_DIRS "${LAHUTA_BOOST_VENDOR_DIR}")
else()
  if(DEFINED Boost_INCLUDE_DIRS)
    list(APPEND _LAHUTA_BGL_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})
  endif()
endif()
if(_LAHUTA_BGL_INCLUDE_DIRS)
  set(CMAKE_REQUIRED_INCLUDES "${_LAHUTA_BGL_INCLUDE_DIRS}")
endif()
check_cxx_source_compiles("#include <boost/graph/adjacency_list.hpp>\nint main(){return 0;}" HAVE_BGL_HEADERS)
unset(CMAKE_REQUIRED_INCLUDES)
if(NOT HAVE_BGL_HEADERS)
  message(FATAL_ERROR "Boost Graph headers not found or not visible to the compiler.\n"
                      "Make sure the submodule exists at core/external/boost (git submodule update --init --recursive).\n"
                      "Otherwise, install Boost >= 1.81 and configure with -DLAHUTA_VENDOR_BOOST=OFF.")
endif()

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

# Set install name for foldseek on macOS
if(APPLE)
  set_target_properties(foldseek PROPERTIES INSTALL_NAME_DIR "@rpath")
endif()
add_subdirectory(external/distopia)

#
# Disable Highway SVE targets in distopia by default to avoid SVE compilation
# issues, such as sizeless types as members. SVE can be enabled explicitly with
# -DLAHUTA_DISTOPIA_ENABLE_SVE=ON
#
if(TARGET libdistopia AND NOT LAHUTA_DISTOPIA_ENABLE_SVE)
  set(LAHUTA_HWY_DISABLE_SVE "(HWY_SVE|HWY_SVE2|HWY_SVE_256|HWY_SVE2_128)")
  target_compile_definitions(libdistopia PRIVATE HWY_DISABLED_TARGETS=${LAHUTA_HWY_DISABLE_SVE})
  message(STATUS "distopia: SVE targets disabled")
endif()

add_library(lmdb STATIC
  ${CMAKE_CURRENT_SOURCE_DIR}/external/lmdb/mdb.c
  ${CMAKE_CURRENT_SOURCE_DIR}/external/lmdb/midl.c
)
set_target_properties(lmdb PROPERTIES LINKER_LANGUAGE C)
target_include_directories(lmdb PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/external)

add_library(lahuta_core SHARED
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
  ${CMAKE_CURRENT_SOURCE_DIR}/external/ctpl
  ${CMAKE_CURRENT_SOURCE_DIR}/external/simde
  ${CMAKE_CURRENT_SOURCE_DIR}/external/gemmi/include
  ${CMAKE_CURRENT_SOURCE_DIR}/external/gemmi/third_party
  ${CMAKE_CURRENT_SOURCE_DIR}/external/gemmi/third_party/zlib
  ${CMAKE_CURRENT_SOURCE_DIR}/external/distopia/include
  ${CMAKE_CURRENT_SOURCE_DIR}/external/lmdb
)

target_link_libraries(lahuta_core PUBLIC
  gemmi_cpp
  rdkit
  libdistopia
  lmdb
  ZLIB::ZLIB
  spdlog::spdlog_header_only
  Boost::headers
)

# Mapping library
add_library(lahuta_mapping STATIC
  src/mapping/bitset_mapping.cpp
  src/mapping/mapping_matrix.cpp
  src/mapping/residue_set.cpp
)
target_compile_features(lahuta_mapping PUBLIC cxx_std_17)
set_property(TARGET lahuta_mapping PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(lahuta_mapping PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/external
)
# SIMDe headers for mapping and dependent code
target_include_directories(lahuta_mapping PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/external/simde
)
# Foldseek headers needed by mapping
target_include_directories(lahuta_mapping PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/external/foldseek
  ${CMAKE_CURRENT_SOURCE_DIR}/external/foldseek/simd
  ${CMAKE_CURRENT_SOURCE_DIR}/external/foldseek/generated
)
target_link_libraries(lahuta_mapping PUBLIC
  lahuta_core
  foldseek
  ZLIB::ZLIB
  spdlog::spdlog_header_only
  Boost::headers
)

# Eigen, system first, bundled fallback
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

if(LAHUTA_BUILD_PYTHON)
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)

  # When building outside scikit-build with no prefix default install
  # prefix to the interpreter's site-packages so cmake --install works.
  if(NOT DEFINED SKBUILD AND CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(_PY_PLATLIB "${Python3_SITEARCH}")
    if(NOT _PY_PLATLIB)
      set(_PY_PLATLIB "${Python3_SITELIB}")
    endif()
    if(_PY_PLATLIB)
      file(TO_CMAKE_PATH "${_PY_PLATLIB}" _PY_PLATLIB)
      set(CMAKE_INSTALL_PREFIX "${_PY_PLATLIB}" CACHE PATH "Install prefix" FORCE)
      message(STATUS "Install prefix not specified. Defaulting to Python site-packages: ${CMAKE_INSTALL_PREFIX}")
    endif()
  endif()

  find_package(pybind11 REQUIRED CONFIG)
  if(pybind11_FOUND)
    message(STATUS "pybind11 found. Building Python bindings.")

    set(PY_SRC_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python/src)
    set(PY_PKG_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../interop/python)
    set(PY_PKG_LIB  ${PY_PKG_ROOT}/lahuta/lib)

    # Staging directory for stub generation
    set(LAHUTA_PY_STAGING_DIR "${CMAKE_BINARY_DIR}/wheel/platlib")

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

    set_target_properties(lahuta_ext PROPERTIES
      BUILD_RPATH   "$<$<PLATFORM_ID:Darwin>:@loader_path;@loader_path/..>$<$<PLATFORM_ID:Linux>:$ORIGIN;$ORIGIN/..>"
      INSTALL_RPATH "$<$<PLATFORM_ID:Darwin>:@loader_path>$<$<PLATFORM_ID:Linux>:$ORIGIN>"
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
    # Link mapping extension against the mapping library so we pull foldseek (which we should only need here)
    target_link_libraries(mapping_ext PRIVATE lahuta_mapping spdlog::spdlog_header_only ZLIB::ZLIB)
    target_compile_features(mapping_ext PRIVATE cxx_std_17)
    set_target_properties(mapping_ext PROPERTIES OUTPUT_NAME "mapping")
    set_target_properties(mapping_ext PROPERTIES
      BUILD_RPATH   "$<$<PLATFORM_ID:Darwin>:@loader_path;@loader_path/..>$<$<PLATFORM_ID:Linux>:$ORIGIN;$ORIGIN/..>"
      INSTALL_RPATH "$<$<PLATFORM_ID:Darwin>:@loader_path>$<$<PLATFORM_ID:Linux>:$ORIGIN>"
    )

    # Install compiled extensions and runtime libraries into the Python package tree
    # These are relative to the CMAKE_INSTALL_PREFIX (which could be set to sitepackages above)
    # Note: When using direct CMake install (cmake --install), both compiled artifacts and Python files are installed by CMake.
    # When using scikit-build-core, only compiled artifacts are installed by CMake, Python files are handled by the build backend.
    install(TARGETS lahuta_ext
      LIBRARY DESTINATION "lahuta/lib"
      COMPONENT python
    )

    install(TARGETS mapping_ext
      LIBRARY DESTINATION "lahuta/lib"
      COMPONENT python
    )

    # runtime shared libraries needed by Python extensions
    install(TARGETS lahuta_core foldseek
      LIBRARY DESTINATION "lahuta/lib"
      COMPONENT python
    )

    # Python file installation depends on build context
    if(NOT DEFINED SKBUILD)
      # Direct CMake install: we need to install Python files ourselves
      install(DIRECTORY ${PY_PKG_ROOT}/lahuta/
        DESTINATION "lahuta"
        COMPONENT python
        FILES_MATCHING
          PATTERN "*.py"
          PATTERN "py.typed"
        PATTERN "__pycache__" EXCLUDE
      )
      message(STATUS "Direct CMake install: Installing Python files via CMake")
    else()
      message(STATUS "scikit-build-core install: Python files handled by build backend")
    endif()

    if(LAHUTA_GENERATE_PY_STUBS)
      set(STUB_GEN_SCRIPT "${PY_PKG_ROOT}/tools/generate_stubs_post_install.py")
      if(NOT EXISTS "${STUB_GEN_SCRIPT}")
        message(FATAL_ERROR "Post-install stub generator not found at ${STUB_GEN_SCRIPT}.")
      endif()

      install(CODE "
        message(STATUS \"Generating Python type stubs...\")
        execute_process(
          COMMAND \"${Python3_EXECUTABLE}\" \"${STUB_GEN_SCRIPT}\" \"${CMAKE_INSTALL_PREFIX}\"
          RESULT_VARIABLE STUB_RESULT
          OUTPUT_VARIABLE STUB_OUTPUT
          ERROR_VARIABLE  STUB_ERROR
          TIMEOUT 60
        )
        if(NOT STUB_RESULT EQUAL 0)
          message(STATUS \"Stub generation output: \${STUB_OUTPUT}\")
          message(WARNING \"Python type stub generation failed with exit code \${STUB_RESULT}: \${STUB_ERROR}\")
          message(STATUS \"Package installation will continue without type stubs\")
        else()
          message(STATUS \"Python type stubs generated successfully\")
        endif()
      "
      COMPONENT python
      )
    endif()
  else()
    message(STATUS "pybind11 not found. Skipping Python bindings build.")
  endif()
endif()

# Tests
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
