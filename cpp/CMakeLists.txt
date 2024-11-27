cmake_minimum_required(VERSION 3.15...3.22)
project(lahuta VERSION 1.0.0 LANGUAGES CXX C)

include(GNUInstallDirs)  # for CMAKE_INSTALL_LIBDIR, etc

if (DEFINED SKBUILD)  # building with scikit-build-core (pip install)
  # we set wheel.install-dir="/data" in pyproject.toml,
  # so here we need to only set paths to /platlib and /scripts
  set(PYTHON_INSTALL_DIR "${SKBUILD_PLATLIB_DIR}")
endif()

option(BUILD_SHARED_LIBS "Build shared library" ON)
option(USE_PYTHON "Build Python bindings" ON)
# option(USE_ZLIB_NG "Use zlib-ng instead of zlib" OFF)
# option(FETCH_ZLIB_NG "Download, build (statically) and use zlib-ng" ON)
# option(INTERNAL_ZLIB "Use subset of zlib distributed with gemmi" OFF)
option(GENERATE_STUBS "Generate Python type stubs" OFF)
option(STANDALONE_PYTHON_MODULE "Avoid linking Python module to liblahuta_cpp DLL" ON)
option(USE_EIGEN3 "Do not use Eigen3" ON)

if (DEFINED Python_EXECUTABLE)
  set(USE_PYTHON ON CACHE BOOL "Build Python bindings" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_BUILD_RPATH_USE_ORIGIN ON)

# set(SANITIZER "thread")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
include(CheckIncludeFile)

if (DEFINED ENV{CXXFLAGS})
  set(USING_ENV_CXXFLAGS ON CACHE BOOL "" FORCE)
endif()

# If CMAKE_BUILD_TYPE is not set, and CXXFLAGS is not set in the environment,
# set it to Release. This is to avoid the default Debug build type.
if (NOT CMAKE_BUILD_TYPE AND NOT USING_ENV_CXXFLAGS)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

if(NOT DEFINED BOOST_ROOT)
  if(DEFINED ENV{CONDA_PREFIX})
    set(BOOST_ROOT $ENV{CONDA_PREFIX})
    else()
     message(FATAL_ERROR "BOOST_ROOT is not defined. Please set it to the path where Boost is installed.")
    endif()
endif()


add_subdirectory(external/gemmi)
add_subdirectory(external/rdkit)
add_subdirectory(external/foldseek)
add_dependencies(foldseek gemmi_cpp)

add_library(lahuta_cpp STATIC
            src/bonds.cpp src/bond_order.cpp src/convert.cpp src/nsgrid.cpp src/neighbors.cpp src/lahuta.cpp src/visitor.cpp src/struct_unit.cpp src/contacts/hydrogen_bonds.cpp 
            src/contacts/geometry.cpp src/nn.cpp src/contacts/charges.cpp 
            src/contacts/hydrophobic.cpp src/contacts/halogen_bonds.cpp
            src/residues.cpp src/contacts/aromaticity.cpp src/aromatics.cpp src/contacts/ionic.cpp
            src/contacts/metalic.cpp src/contacts/cationpi.cpp src/contacts/pistacking.cpp
            src/entities.cpp
            src/contacts/metals.cpp 
            src/contacts/interactions.cpp src/valence_model.cpp src/contacts/groups.cpp
            src/ob/kekulize.cpp src/ob/clean_mol.cpp src/ob/bitvec.cpp
)

set_property(TARGET lahuta_cpp PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(lahuta_cpp PUBLIC gemmi_cpp rdkit foldseek)

find_package(Eigen3 QUIET)
if(USE_EIGEN3)
  if(NOT EIGEN3_FOUND)
    add_subdirectory(external/Eigen)
    add_definitions(-DLAHUTA_USE_EIGEN=1)
  endif()
endif()

target_include_directories(lahuta_cpp PUBLIC
  ${CMAKE_SOURCE_DIR}/include # our headers
  ${CMAKE_SOURCE_DIR}/external # for external headers
  ${CMAKE_SOURCE_DIR}/external/foldseek
  ${CMAKE_SOURCE_DIR}/external/gemmi/include # for gemmi headers
)

add_library(lahuta_shared_headers INTERFACE)
target_include_directories(lahuta_shared_headers INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/include/fseek
  ${CMAKE_CURRENT_SOURCE_DIR}/external/ctpl
)

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

    target_link_libraries(${target} PUBLIC lahuta_cpp)
endfunction()

# add_executable(lex_parallel main.cpp)
# if (LAHUTA_HAS_EIGEN)
#     target_link_libraries(lex_parallel PRIVATE Eigen3::Eigen)
# endif()
# target_link_libraries(lex_parallel PUBLIC lahuta_cpp) 
add_lahuta_executable(lex main.cpp)
add_lahuta_executable(lex_parallel main_threads.cpp)
add_lahuta_executable(lex_fseek main_fseek.cpp)

### Python bindings ###
if (USE_PYTHON)
  message(STATUS "Building Python module.")
  # CMake >=3.18 has subcomponent Development.Module, scikit-build-core also has it
  if (${CMAKE_VERSION} VERSION_LESS 3.18 AND NOT SKBUILD)
    find_package(Python ${PYTHON_VERSION} REQUIRED COMPONENTS Interpreter Development)
  else()
    find_package(Python ${PYTHON_VERSION} REQUIRED COMPONENTS Interpreter Development.Module)
  endif()

  # use pybind11-config (if available) to determine pybind11_DIR
  execute_process(COMMAND pybind11-config --cmakedir OUTPUT_VARIABLE pybind11_DIR)
  string(STRIP "${pybind11_DIR}" pybind11_DIR)
  find_package(pybind11 2.6 CONFIG REQUIRED)
  message(STATUS "Found pybind11 ${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")

  pybind11_add_module(lahuta_py python/at.cpp python/common.cpp python/contacts.cpp python/lahuta.cpp)

  if (0) # FIX: STANDALONE_PYTHON_MODULE branch is broken
    target_sources(lahuta_py PRIVATE 
    $<TARGET_OBJECTS:lahuta_cpp>
    )
    get_target_property(_lahuta_cpp_libs lahuta_cpp LINK_LIBRARIES)
    set_property(TARGET lahuta_py APPEND PROPERTY LINK_LIBRARIES ${_lahuta_cpp_libs})
  else()
    target_link_libraries(lahuta_py PRIVATE lahuta_cpp)
  endif()

  set_property(TARGET lahuta_py PROPERTY OUTPUT_NAME _lahuta)
  if(CMAKE_CXX_FLAGS MATCHES "-Wshadow")
    target_compile_options(lahuta_py PRIVATE "-Wno-shadow")
  endif()

  if (GENERATE_STUBS AND CMAKE_CROSSCOMPILING)
    message(WARNING "Stubs cannot be generated when cross-compiling - skipping.")
  elseif (GENERATE_STUBS)
    set(sep ":")
    if (WIN32)
      set(sep ";")
    endif()
    add_custom_command(TARGET lahuta_py POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E env
         "PYTHONPATH=$<TARGET_FILE_DIR:lahuta_py>${sep}$ENV{PYTHONPATH}"
       ${Python_EXECUTABLE} -m pybind11_stubgen _lahuta
     COMMAND ${CMAKE_COMMAND} -E touch stubs/py.typed
     COMMENT "Generate type stubs with pybind11-stubgen"
     BYPRODUCTS stubs/__init__.pyi stubs/_lahuta.pyi stubs/py.typed
     VERBATIM)
  endif()
else()
  message(STATUS "Skipping Python module. Add -D USE_PYTHON=1 to build it.")
endif()


### Installation ###
# Python module installation:
set(PYTHON_INSTALL_DIR "${CMAKE_BINARY_DIR}/../../lahuta/lib")
message(STATUS "Install directory for Python module: ${PYTHON_INSTALL_DIR}")


if (USE_PYTHON)
  if (DEFINED PYTHON_INSTALL_DIR)
    message(STATUS "Install directory for Python module: ${PYTHON_INSTALL_DIR}")
    set(Python_SITEARCH "${PYTHON_INSTALL_DIR}")
  endif()
  file(TO_CMAKE_PATH "${Python_SITEARCH}" Python_SITEARCH)
  # Using Python_SITEARCH/SITELIB is not good, because they are absolute
  # and don't respect CMAKE_INSTALL_PREFIX.
  # https://discourse.cmake.org/t/findpython3-how-to-specify-local-installation-directory-for-python-module/3580/5
  if (BUILD_SHARED_LIBS AND NOT STANDALONE_PYTHON_MODULE)
    if (DEFINED SKBUILD)
      # When building via pip, Python_SITEARCH is /tmp.../wheel/platlib,
      # different from the actual install path. Best guess:
      set(rpath_rel "../../../${CMAKE_INSTALL_LIBDIR}")
    else()
      file(RELATIVE_PATH rpath_rel
           "${Python_SITEARCH}" "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
    endif()
    set_target_properties(lahuta_py PROPERTIES INSTALL_RPATH "${origin}/${rpath_rel}")
  endif()
  install(TARGETS lahuta_py DESTINATION "${Python_SITEARCH}" COMPONENT py)
  if (GENERATE_STUBS AND NOT CMAKE_CROSSCOMPILING)
    install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/stubs"
            DESTINATION "${Python_SITEARCH}" COMPONENT py)
  endif()
endif()
