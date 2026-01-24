cmake_minimum_required(VERSION 3.20)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(LahutaVersion)
lahuta_read_version_file()

project(lahuta VERSION ${LAHUTA_VERSION_BASE} LANGUAGES CXX C)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  message(FATAL_ERROR "Configure from the repository root (cmake -S . -B build [...])")
endif()

include(CTest) # Provides BUILD_TESTING option
include(GNUInstallDirs)
include(CheckCXXSourceCompiles)

if(APPLE)
  set(CMAKE_MACOSX_RPATH 1)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

  set(CMAKE_INSTALL_RPATH "@loader_path;@loader_path/..") # relative rpath in installation

  set(CMAKE_SKIP_INSTALL_RPATH FALSE)
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
endif()

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
    set(_LAHUTA_OFLAG "-ffast-math")
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


set(RDK_USE_BOOST_SERIALIZATION OFF CACHE BOOL "Disable rdkit boost serialization" FORCE)

set(SPDLOG_FMT_EXTERNAL      OFF CACHE BOOL "Use bundled fmt for spdlog")
set(SPDLOG_BUILD_EXAMPLES    OFF CACHE BOOL "" )
set(SPDLOG_BUILD_TESTS       OFF CACHE BOOL "" )
set(SPDLOG_BUILD_EXTRAS      OFF CACHE BOOL "" )
set(SPDLOG_BUILD_PERF_TESTS  OFF CACHE BOOL "" )
set(SPDLOG_BUILD_HEADER_ONLY ON  CACHE BOOL "" )

# Do not build extras/tests for Highway
set(HWY_ENABLE_TESTS      OFF CACHE BOOL "Disable highway tests"           FORCE)
set(HWY_ENABLE_EXAMPLES   OFF CACHE BOOL "Disable highway examples"        FORCE)
set(HWY_ENABLE_BENCHMARKS OFF CACHE BOOL "Disable highway benchmarks"      FORCE)
set(HWY_ENABLE_CONTRIB    OFF CACHE BOOL "Disable highway contrib targets" FORCE)
set(HWY_ENABLE_INSTALL    OFF CACHE BOOL "Disable highway install targets" FORCE)


find_package(ZLIB REQUIRED)

# Boost headers
if(LAHUTA_VENDOR_BOOST)
  set(LAHUTA_BOOST_VENDOR_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/boost")
  if(NOT EXISTS "${LAHUTA_BOOST_VENDOR_DIR}/boost/version.hpp")
    message(FATAL_ERROR "Vendored Boost headers not found at ${CMAKE_CURRENT_SOURCE_DIR}/external/boost. "
                        "Ensure submodules are initialized (git submodule update --init --recursive), "
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
# Make sure the probe compiles as C++17 regardless of compiler defaults.
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
set(_LAHUTA_SAVE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
if(MSVC)
  string(APPEND CMAKE_REQUIRED_FLAGS " /std:c++17")
else()
  string(APPEND CMAKE_REQUIRED_FLAGS " -std=c++17")
endif()
check_cxx_source_compiles("#include <boost/graph/adjacency_list.hpp>\nint main(){return 0;}" HAVE_BGL_HEADERS)
set(CMAKE_REQUIRED_FLAGS "${_LAHUTA_SAVE_REQUIRED_FLAGS}")
unset(CMAKE_REQUIRED_INCLUDES)
if(NOT HAVE_BGL_HEADERS)
  message(FATAL_ERROR "Boost Graph headers not found or not visible to the compiler. "
                      "Make sure the submodule exists at core/external/boost (git submodule update --init --recursive). "
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

if(TARGET gemmi_cpp)
  target_compile_features(gemmi_cpp PUBLIC cxx_std_17)
endif()
if(TARGET gemmi_headers)
  target_compile_features(gemmi_headers INTERFACE cxx_std_17)
endif()

add_subdirectory(external/rdkit)

if(LAHUTA_BUILD_FOLDSEEK)
  add_subdirectory(external/foldseek)
  add_dependencies(foldseek gemmi_cpp)

  # Set install name for foldseek on macOS
  if(APPLE)
    set_target_properties(foldseek PROPERTIES INSTALL_NAME_DIR "@rpath")
  endif()
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

set(_lahuta_core_library_type SHARED)
if(NOT LAHUTA_BUILD_SHARED_CORE)
  set(_lahuta_core_library_type STATIC)
endif()

add_library(lahuta_core ${_lahuta_core_library_type}
  src/bonds/bonds.cpp
  src/bonds/bond_order.cpp
  src/bonds/clean_bonds.cpp
  src/io/convert.cpp
  src/spatial/kd_tree.cpp
  src/spatial/kd_index.cpp
  src/lahuta.cpp
  src/selections/visitor.cpp
  src/utils/struct_unit.cpp
  src/topology.cpp
  src/residues/residues.cpp
  src/rings/aromatics.cpp
  src/properties/luni_props.cpp
  src/logging/logging.cpp
  src/md/element_utils.cpp
  src/md/gro_reader.cpp
  src/md/xdr.cpp
  src/md/xtc_reader.cpp

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
  src/chemistry/valence_model.cpp

  src/contacts/aromaticity.cpp

  src/chemistry/types/charges.cpp
  src/chemistry/types/metalics.cpp
  src/chemistry/types/hydrophobic.cpp
  src/chemistry/types/halogens.cpp
  src/chemistry/types/hbonding.cpp
  src/typing/getcontacts/atom_typing.cpp

  src/topology/engine.cpp
  src/topology/kernels/neighbors.cpp
  src/topology/kernels/bonds.cpp
  src/topology/kernels/rings.cpp
  src/topology/kernels/residues.cpp
  src/topology/kernels/atom_typing.cpp
  src/topology/kernels/seed_from_model.cpp

  src/bonds/perception/consistency.cpp
  src/bonds/perception/template.cpp
  src/bonds/perception/subgraph.cpp
  src/bonds/perception/nonstandard_bonds.cpp
  src/bonds/perception/subset_merge.cpp

  src/chemistry/neighbors.cpp
  src/chemistry/geometry.cpp
  src/chemistry/predicates.cpp

  src/spatial/nsresults.cpp
  src/spatial/fastns.cpp
  src/spatial/nsutils.cpp

  src/distances/api.cpp
  src/distances/kernels.cpp
  src/distances/matrix.cpp
  src/distances/neighbors.cpp

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

  src/contacts/getcontacts/hbonds.cpp
  src/contacts/getcontacts/hydrophobic.cpp
  src/contacts/getcontacts/salt_bridges.cpp
  src/contacts/getcontacts/vanderwaals.cpp
  src/contacts/getcontacts/pi_cation.cpp
  src/contacts/getcontacts/pistacking.cpp

)

target_compile_features(lahuta_core PUBLIC cxx_std_17)
set_property(TARGET lahuta_core PROPERTY POSITION_INDEPENDENT_CODE ON)

lahuta_configure_core_version(lahuta_core)

target_include_directories(lahuta_core PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/external
  ${CMAKE_CURRENT_SOURCE_DIR}/external/rdkit
  ${CMAKE_CURRENT_SOURCE_DIR}/external/ctpl
  ${CMAKE_CURRENT_SOURCE_DIR}/external/simde
  ${CMAKE_CURRENT_SOURCE_DIR}/external/gemmi/include
  ${CMAKE_CURRENT_SOURCE_DIR}/external/gemmi/third_party
  ${CMAKE_CURRENT_SOURCE_DIR}/external/gemmi/third_party/zlib
  ${CMAKE_CURRENT_SOURCE_DIR}/external/distopia/include
  ${CMAKE_CURRENT_SOURCE_DIR}/external/lmdb
)

#
# Require GCC 9.1+ for full C++17 filesystem support without separate linking
# Not tested.
#
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.1")
    message(FATAL_ERROR "GCC ${CMAKE_CXX_COMPILER_VERSION} detected. Lahuta requires GCC 9.1 or later.")
  endif()
endif()

set(_lahuta_core_private_libs
  rdkit
  libdistopia
  ZLIB::ZLIB
)
set(_lahuta_core_public_libs
  gemmi_cpp
  lmdb
  spdlog::spdlog_header_only
  Boost::headers
)

if(LAHUTA_BUILD_SHARED_CORE)
  target_link_libraries(lahuta_core PRIVATE ${_lahuta_core_private_libs})
  target_link_libraries(lahuta_core PUBLIC ${_lahuta_core_public_libs})
else()
  target_link_libraries(lahuta_core PUBLIC ${_lahuta_core_private_libs} ${_lahuta_core_public_libs})
endif()

# Mapping library
if(LAHUTA_BUILD_FOLDSEEK)
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
endif()

if(LAHUTA_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(LAHUTA_VENDOR_EIGEN)
  set(LAHUTA_EIGEN_VENDOR_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/Eigen")
  if(NOT EXISTS "${LAHUTA_EIGEN_VENDOR_DIR}/eigen/Eigen/Core")
    message(FATAL_ERROR "Vendored Eigen headers not found at ${CMAKE_CURRENT_SOURCE_DIR}/external/Eigen. "
                        "Ensure submodules are initialized (git submodule update --init --recursive), "
                        "or disable vendoring with -DLAHUTA_VENDOR_EIGEN=OFF to use system Eigen.")
  endif()
  add_subdirectory(external/Eigen)
  message(STATUS "Using vendored Eigen headers at ${LAHUTA_EIGEN_VENDOR_DIR}")
else()
  find_package(Eigen3 REQUIRED)
  if(NOT TARGET Eigen3::Eigen)
    message(FATAL_ERROR "System Eigen3 not found. Install Eigen3 or configure with -DLAHUTA_VENDOR_EIGEN=ON to use vendored Eigen.")
  endif()
  message(STATUS "Using system Eigen3")
endif()

target_link_libraries(lahuta_core PUBLIC Eigen3::Eigen)
target_compile_definitions(lahuta_core PUBLIC LAHUTA_USE_EIGEN=1)

if(LAHUTA_BUILD_PYTHON)
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
  find_package(Python3 COMPONENTS NumPy)

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
      ${PY_SRC_DIR}/elements.cpp
      ${PY_SRC_DIR}/records.cpp
      ${PY_SRC_DIR}/properties.cpp

      ${PY_SRC_DIR}/grid_search.cpp
      ${PY_SRC_DIR}/distance.cpp
      ${PY_SRC_DIR}/contacts.cpp

      ${PY_SRC_DIR}/rdkit.cpp

      ${PY_SRC_DIR}/pipeline_dynamic.cpp
      ${PY_SRC_DIR}/db.cpp
    )

    pybind11_add_module(lahuta_ext MODULE ${LAHUTA_PY_SOURCES})
    target_include_directories(lahuta_ext PRIVATE ${PY_SRC_DIR})
    target_link_libraries(lahuta_ext PRIVATE lahuta_core spdlog::spdlog_header_only ZLIB::ZLIB)
    if(TARGET Python3::NumPy)
      target_link_libraries(lahuta_ext PRIVATE Python3::NumPy)
      message(STATUS "NumPy C API found. Enabling optimized array operations for lahuta_ext.")
    else()
      message(STATUS "NumPy C API not found. Using pybind11 fallback for array operations.")
    endif()
    target_compile_features(lahuta_ext PRIVATE cxx_std_17)
    set_target_properties(lahuta_ext PROPERTIES OUTPUT_NAME "lahuta")

    set_target_properties(lahuta_ext PROPERTIES
      BUILD_RPATH   "$<$<PLATFORM_ID:Darwin>:@loader_path;@loader_path/..>$<$<PLATFORM_ID:Linux>:$ORIGIN;$ORIGIN/..>"
      INSTALL_RPATH "$<$<PLATFORM_ID:Darwin>:@loader_path;@loader_path/..>$<$<PLATFORM_ID:Linux>:$ORIGIN;$ORIGIN/..>"
      # Build with install RPATH to prevent re-install errors on macOS
      BUILD_WITH_INSTALL_RPATH TRUE
    )

    # mapping submodule, depends at Python-import level on lahuta
    if(LAHUTA_BUILD_FOLDSEEK)
      set(LAHUTA_MAPPING_SOURCES
        ${PY_SRC_DIR}/mapping/mapping.cpp
        ${PY_SRC_DIR}/mapping/fseek.cpp
        ${PY_SRC_DIR}/mapping/align.cpp
        ${PY_SRC_DIR}/mapping/te.cpp
      )
      pybind11_add_module(mapping_ext MODULE ${LAHUTA_MAPPING_SOURCES})
      target_include_directories(mapping_ext PRIVATE ${PY_SRC_DIR})
      target_link_libraries(mapping_ext PRIVATE lahuta_mapping spdlog::spdlog_header_only ZLIB::ZLIB)
      if(TARGET Python3::NumPy)
        target_link_libraries(mapping_ext PRIVATE Python3::NumPy)
      endif()
      target_compile_features(mapping_ext PRIVATE cxx_std_17)
      set_target_properties(mapping_ext PROPERTIES OUTPUT_NAME "mapping")
      set_target_properties(mapping_ext PROPERTIES
        BUILD_RPATH   "$<$<PLATFORM_ID:Darwin>:@loader_path;@loader_path/..>$<$<PLATFORM_ID:Linux>:$ORIGIN;$ORIGIN/..>"
        INSTALL_RPATH "$<$<PLATFORM_ID:Darwin>:@loader_path;@loader_path/..>$<$<PLATFORM_ID:Linux>:$ORIGIN;$ORIGIN/..>"
        # Build with install RPATH to prevent re-install errors on macOS
        BUILD_WITH_INSTALL_RPATH TRUE
      )
    endif()

    # Install compiled extensions and runtime libraries into the Python package tree
    # These are relative to the CMAKE_INSTALL_PREFIX (which could be set to sitepackages above)
    # Note: When using direct CMake install (cmake --install), both compiled artifacts and Python files are installed by CMake.
    # When using scikit-build-core, only compiled artifacts are installed by CMake, Python files are handled by the build backend.
    install(TARGETS lahuta_ext
      LIBRARY DESTINATION "lahuta/lib"
      COMPONENT python
    )

    if(LAHUTA_BUILD_FOLDSEEK)
      install(TARGETS mapping_ext
        LIBRARY DESTINATION "lahuta/lib"
        COMPONENT python
      )
    endif()

    # runtime shared libraries needed by Python extensions
    if(LAHUTA_BUILD_SHARED_CORE)
      set(_lahuta_core_python_dest LIBRARY DESTINATION "lahuta/lib")
    else()
      set(_lahuta_core_python_dest ARCHIVE DESTINATION "lahuta/lib")
    endif()

    if(LAHUTA_BUILD_FOLDSEEK)
      install(TARGETS lahuta_core
        ${_lahuta_core_python_dest}
        COMPONENT python
      )
      install(TARGETS foldseek
        LIBRARY DESTINATION "lahuta/lib"
        COMPONENT python
      )
    else()
      install(TARGETS lahuta_core
        ${_lahuta_core_python_dest}
        COMPONENT python
      )
    endif()

    # Configure Python version file for direct CMake installs
    lahuta_configure_python_version(
      "${PY_PKG_ROOT}"
      "${CMAKE_CURRENT_BINARY_DIR}/python_configured"
    )

    # Python file installation depends on build context
    if(NOT DEFINED SKBUILD)
      # Direct CMake install: we need to install Python files ourselves
      install(DIRECTORY ${PY_PKG_ROOT}/lahuta/
        DESTINATION "lahuta"
        COMPONENT python
        FILES_MATCHING
          PATTERN "*.py"
          PATTERN "*.pyi"
          PATTERN "py.typed"
        PATTERN "__pycache__" EXCLUDE
        PATTERN "_version.py" EXCLUDE  # Exclude template, use configured version below
      )
      # Install the CMake-configured version file
      install(FILES "${LAHUTA_PYTHON_VERSION_FILE}"
        DESTINATION "lahuta"
        COMPONENT python
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
