cmake_minimum_required(VERSION 3.15...3.22)
project(lahuta VERSION 1.0.0 LANGUAGES CXX C)

include(GNUInstallDirs)  # for CMAKE_INSTALL_LIBDIR, etc

if (DEFINED SKBUILD)  # building with scikit-build-core (pip install)
  # we set wheel.install-dir="/data" in pyproject.toml,
  # so here we need to only set paths to /platlib and /scripts
  set(PYTHON_INSTALL_DIR "${SKBUILD_PLATLIB_DIR}")
endif()

option(USE_PYTHON "Build Python bindings" ON)
option(USE_ZLIB_NG "Use zlib-ng instead of zlib" ON)
option(FETCH_ZLIB_NG "Download, build (statically) and use zlib-ng" OFF)
option(INTERNAL_ZLIB "Use subset of zlib distributed with gemmi" OFF)
option(GENERATE_STUBS "Generate Python type stubs" OFF)
option(STANDALONE_PYTHON_MODULE "Avoid linking Python module to liblahuta_cpp DLL" ON)

if (DEFINED Python_EXECUTABLE)
  set(USE_PYTHON ON CACHE BOOL "Build Python bindings" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_BUILD_RPATH_USE_ORIGIN ON)

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

if(NOT DEFINED RDKit_ROOT)
    if(DEFINED ENV{CONDA_PREFIX})
        set(RDKit_ROOT $ENV{CONDA_PREFIX})
    else()
        message(FATAL_ERROR "No rkdit")
    endif()
endif()

list(APPEND CMAKE_PREFIX_PATH ${RDKit_ROOT}) 

# Find external dependencies
find_package(RDKit REQUIRED)


### Compile gemmi in external/gemmi ###
add_subdirectory(external/gemmi)




file(GLOB_RECURSE C_SOURCES
    "external/gemmi/third_party/zlib/*.c"
)
#
add_library(lahuta_cpp SHARED
            src/bonds.cpp src/bond_order.cpp src/conv.cpp src/nsgrid.cpp 
            src/bond_table/bonds.cpp src/ob/kekulize.cpp src/ob/bond_utils.cpp 
            # external/gemmi/src/ddl.cpp external/gemmi/src/gz.cpp external/gemmi/src/json.cpp
            # external/gemmi/src/mmcif.cpp external/gemmi/src/mmread_gz.cpp external/gemmi/src/polyheur.cpp
            # external/gemmi/src/read_cif.cpp external/gemmi/src/resinfo.cpp
            external/ob/bitvec.cpp external/ob/elements.cpp 
            # ${C_SOURCES}
)

set_property(TARGET lahuta_cpp PROPERTY POSITION_INDEPENDENT_CODE ON)
# target_compile_definitions(lahuta_cpp PRIVATE LAHUTA_BUILD)

target_include_directories(lahuta_cpp PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/external/ob
  # ${CMAKE_CURRENT_SOURCE_DIR}/external/gemmi/include
  # ${CMAKE_CURRENT_SOURCE_DIR}/external/gemmi/third_party
  # ${CMAKE_CURRENT_SOURCE_DIR}/external/gemmi/third_party/zlib
)

target_link_libraries(lahuta_cpp PUBLIC
    $ENV{CONDA_PREFIX}/lib/libRDKitGraphMol.so
    $ENV{CONDA_PREFIX}/lib/libRDKitRDGeneral.so
    $ENV{CONDA_PREFIX}/lib/libRDKitRDGeometryLib.so
    $ENV{CONDA_PREFIX}/lib/libRDKitSubstructMatch.so
    $ENV{CONDA_PREFIX}/lib/libRDKitSmilesParse.so
)


target_compile_definitions(lahuta_cpp PUBLIC GEMMI_SHARED)

# target_sources(lahuta_cpp PRIVATE
#       "external/gemmi/third_party/zlib/adler32.c"
#       "external/gemmi/third_party/zlib/crc32.c"
#       "external/gemmi/third_party/zlib/gzlib.c"
#       "external/gemmi/third_party/zlib/gzread.c"
#       "external/gemmi/third_party/zlib/inflate.c"
#       "external/gemmi/third_party/zlib/inftrees.c"
#       "external/gemmi/third_party/zlib/inffast.c"
#       "external/gemmi/third_party/zlib/zutil.c")
# check_include_file(unistd.h has_unistd_h)
# target_compile_definitions(lahuta_cpp PRIVATE NO_GZCOMPRESS=1 DYNAMIC_CRC_TABLE=1)
# if (has_unistd_h)
#   target_compile_definitions(lahuta_cpp PRIVATE Z_HAVE_UNISTD_H=1)
# endif()
# target_include_directories(lahuta_cpp PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/external/gemmi/third_party/zlib")
# if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
#   # /wd4267 - conversion from 'size_t' to 'unsigned int', possible loss of data
#   target_compile_options(lahuta_cpp PRIVATE /wd4267)
# endif()
#
# link lahuta_cpp to gemmi
# target_link_libraries(lahuta_cpp PRIVATE gemmi_cpp)
# target_link_libraries(lahuta_cpp PRIVATE gemmi_cpp)

# merge gemmi and lahuta_cpp into one library
# add_library(lahuta::lahuta_cpp ALIAS lahuta_cpp)

# Define the library target
# add_library(lahuta_cpp SHARED ${LAHUTA_SOURCES})
target_link_libraries(lahuta_cpp PRIVATE gemmi_cpp)
#
# target_include_directories(lahuta_cpp PUBLIC
#   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#   $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
# )

target_compile_definitions(lahuta_cpp PRIVATE LAHUTA_BUILD)
# set_target_properties(lahuta_cpp PROPERTIES
#   VERSION ${PROJECT_VERSION}
#   SOVERSION ${PROJECT_VERSION_MAJOR}
#   OUTPUT_NAME lahuta_cpp
# )
#
# # Install the library
# install(TARGETS lahuta_cpp
#   EXPORT lahutaTargets
#   LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
#   ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
# )
#
# # Install the headers
# # install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
#
# # Export the targets
# install(EXPORT lahutaTargets
#   FILE lahutaTargets.cmake
#   NAMESPACE lahuta::
#   DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/lahuta
# )





# Define the executable target
add_executable(lex main.cpp)
# target_link_libraries(lex PRIVATE lahuta_cpp)
# link to gemmi instead
target_link_libraries(lex PRIVATE lahuta_cpp gemmi_cpp)


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

  pybind11_add_module(lahuta_py python/lahuta.cpp)
  target_link_libraries(lahuta_py PRIVATE lahuta_cpp gemmi_cpp)

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
       ${Python_EXECUTABLE} -m pybind11_stubgen gemmi
         --enum-class-locations=Ignore:gemmi.ContactSearch
         --enum-class-locations=.+:gemmi
         --numpy-array-remove-parameters
         --output-dir stubs
     COMMAND ${CMAKE_COMMAND} -E touch stubs/gemmi/py.typed
     COMMENT "Generate type stubs with pybind11-stubgen"
     BYPRODUCTS stubs/gemmi/__init__.pyi stubs/gemmi/cif.pyi stubs/gemmi/py.typed
     VERBATIM)
  endif()
else()
  message(STATUS "Skipping Python module. Add -D USE_PYTHON=1 to build it.")
endif()


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
  install(DIRECTORY examples DESTINATION "${Python_SITEARCH}/lahuta" COMPONENT py
          FILES_MATCHING PATTERN "*.py"
          PATTERN "[._]*" EXCLUDE)
  if (GENERATE_STUBS AND NOT CMAKE_CROSSCOMPILING)
    install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/stubs/lahuta"
            DESTINATION "${Python_SITEARCH}" COMPONENT py)
  endif()
endif()

