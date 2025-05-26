cmake_minimum_required(VERSION 3.15...3.22)
project(lahuta VERSION 1.0.0 LANGUAGES CXX C)

include(GNUInstallDirs)  # for CMAKE_INSTALL_LIBDIR, etc

set(HWY_ENABLE_TESTS OFF CACHE BOOL "Enable HWY tests" FORCE)
set(HWY_ENABLE_EXAMPLES OFF CACHE BOOL "Disable highway examples" FORCE)
set(HWY_ENABLE_BENCHMARKS OFF CACHE BOOL "Disable highway benchmarks" FORCE)
set(HWY_ENABLE_CONTRIB OFF CACHE BOOL "Disable highway contrib targets" FORCE)
set(HWY_ENABLE_INSTALL OFF CACHE BOOL "Disable highway install targets" FORCE)

set(SPDLOG_FMT_EXTERNAL OFF CACHE BOOL "Use bundled fmt")
set(SPDLOG_BUILD_EXAMPLES OFF  CACHE BOOL "")
set(SPDLOG_BUILD_TESTS OFF     CACHE BOOL "")
set(SPDLOG_BUILD_EXTRAS OFF    CACHE BOOL "")
set(SPDLOG_BUILD_PERF_TESTS OFF CACHE BOOL "")
set(SPDLOG_BUILD_HEADER_ONLY ON CACHE BOOL "")

if (DEFINED SKBUILD)  # building with scikit-build-core (pip install)
  # we set wheel.install-dir="/data" in pyproject.toml,
  # so here we need to only set paths to /platlib and /scripts
  set(PYTHON_INSTALL_DIR "${SKBUILD_PLATLIB_DIR}")
endif()

option(BUILD_SHARED_LIBS "Build shared library" ON)
option(USE_PYTHON "Build Python bindings" ON)

option(USE_ZLIB_NG "Use zlib-ng instead of zlib" OFF)
option(FETCH_ZLIB_NG "Download, build (statically) and use zlib-ng" OFF)
option(INTERNAL_ZLIB "Use subset of zlib distributed with gemmi" OFF)
option(GENERATE_STUBS "Generate Python type stubs" OFF)
option(STANDALONE_PYTHON_MODULE "Avoid linking Python module to liblahuta_cpp DLL" ON)
option(USE_EIGEN3 "Do not use Eigen3" ON)

if (DEFINED Python_EXECUTABLE)
  set(USE_PYTHON ON CACHE BOOL "Build Python bindings" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_BUILD_RPATH_USE_ORIGIN ON)

# set(SANITIZER "thread")

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release)
endif()

# set(CMAKE_BUILD_TYPE Release)
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native -flto")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -flto -ffast-math")
# set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# find_library(JEMALLOC_LIB jemalloc REQUIRED)

include(CheckIncludeFile)

if (DEFINED ENV{CXXFLAGS})
  set(USING_ENV_CXXFLAGS ON CACHE BOOL "" FORCE)
endif()

# If CMAKE_BUILD_TYPE is not set, and CXXFLAGS is not set in the environment,
# set it to Release. This is to avoid the default Debug build type.
#if (NOT CMAKE_BUILD_TYPE AND NOT USING_ENV_CXXFLAGS)
#  set(CMAKE_BUILD_TYPE Release CACHE STRING
#      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
#      FORCE)
#endif()
#message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# We need to find zlib before adding gemmi so it uses the same zlib
find_package(ZLIB REQUIRED)
if(ZLIB_FOUND)
  message(STATUS "Found zlib: ${ZLIB_LIBRARIES}")
else()
  message(FATAL_ERROR "zlib not found - required for compression")
endif()

find_package(Boost REQUIRED COMPONENTS iostreams)

if(Boost_FOUND)
  message(STATUS "Boost found at ${Boost_INCLUDE_DIRS}")
else()
  message(FATAL_ERROR "Boost not found!")
endif()


if(NOT DEFINED BOOST_ROOT)
  if(DEFINED ENV{CONDA_PREFIX})
    set(BOOST_ROOT $ENV{CONDA_PREFIX})
    else()
     message(FATAL_ERROR "BOOST_ROOT is not defined. Please set it to the path where Boost is installed.")
    endif()
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

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

add_library(lahuta_cpp STATIC
            src/bonds.cpp
            src/bond_order.cpp
            src/convert.cpp
            src/nsgrid.cpp
            src/lahuta.cpp
            src/visitor.cpp
            src/struct_unit.cpp
            src/contacts/hydrogen_bonds.cpp
            src/contacts/charges.cpp
            src/topology.cpp
            src/contacts/hydrophobic.cpp
            src/contacts/halogen_bonds.cpp
            src/residues.cpp
            src/contacts/aromaticity.cpp
            src/aromatics.cpp
            src/luni_props.cpp
            src/logging.cpp
            src/models/topology.cpp
            src/contacts/metals.cpp
            src/models/parser.cpp
            external/mmap/MemoryMapped.cpp
            src/valence_model.cpp
            src/contacts/groups.cpp
            src/ob/kekulize.cpp
            src/ob/clean_mol.cpp
            src/ob/bitvec.cpp

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

)

set_property(TARGET lahuta_cpp PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(lahuta_cpp PUBLIC gemmi_cpp rdkit foldseek libdistopia lmdb ZLIB::ZLIB spdlog::spdlog_header_only)

target_include_directories(lahuta_cpp PUBLIC
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/external
  ${CMAKE_SOURCE_DIR}/external/foldseek
  ${CMAKE_SOURCE_DIR}/external/gemmi/include
  ${CMAKE_SOURCE_DIR}/external/gemmi/third_party   # For zlib headers
  ${CMAKE_SOURCE_DIR}/external/gemmi/third_party/zlib
  ${CMAKE_SOURCE_DIR}/external/distopia/include
  ${CMAKE_SOURCE_DIR}/external/lmdb
)

add_library(lahuta_shared_headers INTERFACE)
target_include_directories(lahuta_shared_headers INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/include/fseek
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

# add_compile_options(-fsanitize=${SANITIZER} -fno-omit-frame-pointer)
# add_link_options(-fsanitize=${SANITIZER} -fno-omit-frame-pointer)

function(add_lahuta_executable target main_file)
    add_executable(${target} ${main_file})

    if (LAHUTA_HAS_EIGEN)
        target_link_libraries(${target} PRIVATE Eigen3::Eigen)
    endif()

    target_link_libraries(${target} PUBLIC lahuta_cpp ZLIB::ZLIB spdlog::spdlog_header_only)
endfunction()

# add_executable(lex_parallel main.cpp)
# if (LAHUTA_HAS_EIGEN)
#     target_link_libraries(lex_parallel PRIVATE Eigen3::Eigen)
# endif()
# target_link_libraries(lex_parallel PUBLIC lahuta_cpp) 


add_lahuta_executable(lex main.cpp)
add_lahuta_executable(lex_ag main_ag.cpp)
add_lahuta_executable(lex_parallel main_threads.cpp)
# add_lahuta_executable(lex_fseek main_fseek.cpp)
# add_lahuta_executable(lex_fut main_future.cpp)
# add_lahuta_executable(lex_align main_align.cpp)
# add_lahuta_executable(create_db create_db.cpp)
# add_lahuta_executable(read_db read_db.cpp)

# ### Python bindings ###
# if (USE_PYTHON)
#   message(STATUS "Building Python module.")
#   # CMake >=3.18 has subcomponent Development.Module, scikit-build-core also has it
#   if (${CMAKE_VERSION} VERSION_LESS 3.18 AND NOT SKBUILD)
#     find_package(Python ${PYTHON_VERSION} REQUIRED COMPONENTS Interpreter Development)
#   else()
#     find_package(Python ${PYTHON_VERSION} REQUIRED COMPONENTS Interpreter Development.Module)
#   endif()
# 
#   # use pybind11-config (if available) to determine pybind11_DIR
#   execute_process(COMMAND pybind11-config --cmakedir OUTPUT_VARIABLE pybind11_DIR)
#   string(STRIP "${pybind11_DIR}" pybind11_DIR)
#   find_package(pybind11 2.6 CONFIG REQUIRED)
#   message(STATUS "Found pybind11 ${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
# 
#   pybind11_add_module(lahuta_py python/at.cpp python/array.cpp python/logger_py.cpp python/contacts.cpp python/py_properties.cpp python/py_align.cpp python/py_te.cpp python/lahuta.cpp)
# 
#   if (0) # FIX: STANDALONE_PYTHON_MODULE branch is broken
#     target_sources(lahuta_py PRIVATE 
#     $<TARGET_OBJECTS:lahuta_cpp>
#     )
#     get_target_property(_lahuta_cpp_libs lahuta_cpp LINK_LIBRARIES)
#     set_property(TARGET lahuta_py APPEND PROPERTY LINK_LIBRARIES ${_lahuta_cpp_libs})
#   else()
#     target_link_libraries(lahuta_py PRIVATE lahuta_cpp)
#   endif()
# 
#   set_property(TARGET lahuta_py PROPERTY OUTPUT_NAME _lahuta)
#   if(CMAKE_CXX_FLAGS MATCHES "-Wshadow")
#     target_compile_options(lahuta_py PRIVATE "-Wno-shadow")
#   endif()
# 
#   if (GENERATE_STUBS AND CMAKE_CROSSCOMPILING)
#     message(WARNING "Stubs cannot be generated when cross-compiling - skipping.")
#   elseif (GENERATE_STUBS)
#     set(sep ":")
#     if (WIN32)
#       set(sep ";")
#     endif()
#     add_custom_command(TARGET lahuta_py POST_BUILD
#     COMMAND ${CMAKE_COMMAND} -E env
#          "PYTHONPATH=$<TARGET_FILE_DIR:lahuta_py>${sep}$ENV{PYTHONPATH}"
#        ${Python_EXECUTABLE} -m pybind11_stubgen _lahuta
#      COMMAND ${CMAKE_COMMAND} -E touch stubs/py.typed
#      COMMENT "Generate type stubs with pybind11-stubgen"
#      BYPRODUCTS stubs/__init__.pyi stubs/_lahuta.pyi stubs/py.typed
#      VERBATIM)
#   endif()
# else()
#   message(STATUS "Skipping Python module. Add -D USE_PYTHON=1 to build it.")
# endif()
# 
# 
# ### Installation ###
# # Python module installation:
# set(PYTHON_INSTALL_DIR "${CMAKE_BINARY_DIR}/../../lahuta/lib")
# message(STATUS "Install directory for Python module: ${PYTHON_INSTALL_DIR}")
# 
# 
# if (USE_PYTHON)
#   if (DEFINED PYTHON_INSTALL_DIR)
#     message(STATUS "Install directory for Python module: ${PYTHON_INSTALL_DIR}")
#     set(Python_SITEARCH "${PYTHON_INSTALL_DIR}")
#   endif()
#   file(TO_CMAKE_PATH "${Python_SITEARCH}" Python_SITEARCH)
#   # Using Python_SITEARCH/SITELIB is not good, because they are absolute
#   # and don't respect CMAKE_INSTALL_PREFIX.
#   # https://discourse.cmake.org/t/findpython3-how-to-specify-local-installation-directory-for-python-module/3580/5
#   if (BUILD_SHARED_LIBS AND NOT STANDALONE_PYTHON_MODULE)
#     if (DEFINED SKBUILD)
#       # When building via pip, Python_SITEARCH is /tmp.../wheel/platlib,
#       # different from the actual install path. Best guess:
#       set(rpath_rel "../../../${CMAKE_INSTALL_LIBDIR}")
#     else()
#       file(RELATIVE_PATH rpath_rel "${Python_SITEARCH}" "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
#     endif()
#     # set_target_properties(lahuta_py PROPERTIES INSTALL_RPATH "${origin}/${rpath_rel}")
#     set_target_properties(lahuta_py PROPERTIES INSTALL_RPATH "@loader_path" BUILD_WITH_INSTALL_RPATH TRUE)
#   endif()
#   install(TARGETS lahuta_py DESTINATION "${Python_SITEARCH}" COMPONENT py)
#   install(FILES "${CMAKE_BINARY_DIR}/external/foldseek/libfoldseek.dylib" DESTINATION "${Python_SITEARCH}" COMPONENT py)
#   install(FILES "${CMAKE_BINARY_DIR}/external/distopia/libdistopia.dylib" DESTINATION "${Python_SITEARCH}" COMPONENT py)
#   if (GENERATE_STUBS AND NOT CMAKE_CROSSCOMPILING)
#     install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/stubs" DESTINATION "${Python_SITEARCH}" COMPONENT py)
#   endif()
# endif()
# 
# add_custom_command(TARGET lahuta_py POST_BUILD
#   COMMAND install_name_tool -add_rpath "@loader_path" $<TARGET_FILE:lahuta_py>
#   COMMENT "Adding @loader_path rpath to _lahuta module"
# )
