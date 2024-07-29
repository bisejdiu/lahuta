cmake_minimum_required(VERSION 3.15...3.22)
project(lahuta VERSION 1.0.0 LANGUAGES CXX C)

include(GNUInstallDirs)  # for CMAKE_INSTALL_LIBDIR, etc

if (DEFINED SKBUILD)  # building with scikit-build-core (pip install)
  # we set wheel.install-dir="/data" in pyproject.toml,
  # so here we need to only set paths to /platlib and /scripts
  set(PYTHON_INSTALL_DIR "${SKBUILD_PLATLIB_DIR}")
  #set(CMAKE_INSTALL_BINDIR "${SKBUILD_SCRIPTS_DIR}")
endif()

# option(BUILD_GEMMI_PROGRAM "Build gemmi command-line program" ON)
option(USE_PYTHON "Build Python bindings" ON)
option(USE_ZLIB_NG "Use zlib-ng instead of zlib" OFF)
option(FETCH_ZLIB_NG "Download, build (statically) and use zlib-ng" OFF)
option(INTERNAL_ZLIB "Use subset of zlib distributed with gemmi" ON)
option(GENERATE_STUBS "Generate Python type stubs" OFF)
option(EXTRA_WARNINGS "Set extra warning flags" OFF)
# option(USE_WMAIN "(Windows only) take Unicode arguments in gemmi program" ON)
option(STANDALONE_PYTHON_MODULE "Avoid linking Python module to liblahuta_cpp DLL" ON)

if (DEFINED Python_EXECUTABLE)
  set(USE_PYTHON ON CACHE BOOL "Build Python bindings" FORCE)
endif()

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

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_BUILD_RPATH_USE_ORIGIN ON)

# build type
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
    # set(BUILD_SHARED_LIBS_here ${BUILD_SHARED_LIBS})
    # set(BUILD_SHARED_LIBS OFF CACHE INTERNAL "build static zlib-ng")
    add_subdirectory(${zlib-ng_SOURCE_DIR} ${zlib-ng_BINARY_DIR} EXCLUDE_FROM_ALL)
    # set(BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS_here} CACHE BOOL "Build shared library" FORCE)
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
    # FIX: Test this path
    set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/external/gemmi/third_party")
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

# Set default RDKit_ROOT if not provided
# FIX: installation should also support pip
if(NOT DEFINED RDKit_ROOT)
    if(DEFINED ENV{CONDA_PREFIX})
        set(RDKit_ROOT $ENV{CONDA_PREFIX})
    else()
        message(FATAL_ERROR "RDKit_ROOT not set and CONDA_PREFIX environment variable not found. Please set RDKit_ROOT manually.")
    endif()
endif()

# Add RDKit_ROOT to CMAKE_PREFIX_PATH
# FIX: should we instead use find_package(RDKit) here? 
list(APPEND CMAKE_PREFIX_PATH ${RDKit_ROOT}) 

# Find external dependencies
find_package(Eigen3 REQUIRED)
find_package(RDKit REQUIRED)

add_library(lahuta
            src/bonds.cpp src/bond_order.cpp src/conv.cpp src/nsgrid.cpp
            src/bond_table/bonds.cpp src/ob/kekulize.cpp src/ob/bond_utils.cpp 
            external/gemmi/src/ddl.cpp external/gemmi/src/gz.cpp external/gemmi/src/json.cpp
            external/gemmi/src/mmcif.cpp external/gemmi/src/mmread_gz.cpp external/gemmi/src/polyheur.cpp
            external/gemmi/src/read_cif.cpp external/gemmi/src/resinfo.cpp
            external/ob/bitvec.cpp external/ob/elements.cpp)

target_compile_features(lahuta INTERFACE cxx_std_17)

set_property(TARGET lahuta PROPERTY POSITION_INDEPENDENT_CODE ON)
target_compile_definitions(lahuta PRIVATE LAHUTA_BUILD)

target_include_directories(lahuta PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_include_directories(lahuta PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/external/ob")
target_include_directories(lahuta PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/external/gemmi/include")
target_include_directories(lahuta PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/external/gemmi/third_party")
target_include_directories(lahuta PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/external/gemmi/third_party/zlib")

target_link_libraries(lahuta PRIVATE
    $ENV{CONDA_PREFIX}/lib/libRDKitGraphMol.so
    $ENV{CONDA_PREFIX}/lib/libRDKitRDGeneral.so
    $ENV{CONDA_PREFIX}/lib/libRDKitRDGeometryLib.so
    $ENV{CONDA_PREFIX}/lib/libRDKitSubstructMatch.so
    $ENV{CONDA_PREFIX}/lib/libRDKitSmilesParse.so
)
set_target_properties(lahuta PROPERTIES
    INSTALL_RPATH "${RDKit_ROOT}/lib"
    BUILD_WITH_INSTALL_RPATH TRUE
)

add_executable(lxahuta_cpp main.cpp)
target_link_libraries(lxahuta_cpp PRIVATE lahuta)


# FIX: test the different paths
if (FETCH_ZLIB_NG)
  target_compile_definitions(lahuta PRIVATE USE_ZLIB_NG=1)
  target_link_libraries(lahuta PRIVATE zlib)
elseif (ZLIB_NG_FOUND)
  target_compile_definitions(lahuta PRIVATE USE_ZLIB_NG=1)
  target_link_libraries(lahuta PRIVATE ZLIB_NG::ZLIB_NG)
elseif (ZLIB_FOUND)
  target_link_libraries(lahuta PRIVATE ZLIB::ZLIB)
else()
  target_sources(lahuta PRIVATE
      "external/gemmi/third_party/zlib/adler32.c"
      "external/gemmi/third_party/zlib/crc32.c"
      "external/gemmi/third_party/zlib/gzlib.c"
      "external/gemmi/third_party/zlib/gzread.c"
      "external/gemmi/third_party/zlib/inflate.c"
      "external/gemmi/third_party/zlib/inftrees.c"
      "external/gemmi/third_party/zlib/inffast.c"
      "external/gemmi/third_party/zlib/zutil.c")
  check_include_file(unistd.h has_unistd_h)
  target_compile_definitions(lahuta PRIVATE NO_GZCOMPRESS=1 DYNAMIC_CRC_TABLE=1)
  if (has_unistd_h)
    target_compile_definitions(lahuta PRIVATE Z_HAVE_UNISTD_H=1)
  endif()
  target_include_directories(lahuta PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/external/gemmi/third_party/zlib")
  if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    # /wd4267 - conversion from 'size_t' to 'unsigned int', possible loss of data
    target_compile_options(lahuta PRIVATE /wd4267)
  endif()
endif()

### Python bindings ###
if (USE_PYTHON)
  message(STATUS "The python module will be built.")
  # CMake >=3.18 has subcomponent Development.Module, scikit-build-core also has it
  if (${CMAKE_VERSION} VERSION_LESS 3.18 AND NOT SKBUILD)
    find_package(Python ${PYTHON_VERSION} REQUIRED COMPONENTS Interpreter Development)
  else()
    find_package(Python ${PYTHON_VERSION} REQUIRED COMPONENTS Interpreter Development.Module)
  endif()
  if (EXISTS "${CMAKE_HOME_DIRECTORY}/pybind11")
        message(STATUS "Using ${CMAKE_HOME_DIRECTORY}/pybind11 (internal copy).") # FIX: what's this?
    add_subdirectory(pybind11)
  else()
    # use pybind11-config (if available) to determine pybind11_DIR
    execute_process(COMMAND pybind11-config --cmakedir OUTPUT_VARIABLE pybind11_DIR)
    string(STRIP "${pybind11_DIR}" pybind11_DIR)
    find_package(pybind11 2.6 CONFIG REQUIRED)
    message(STATUS "Found pybind11 ${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
  endif()
  pybind11_add_module(lahuta_py
          python/lahuta.cpp)
  if (STANDALONE_PYTHON_MODULE)
    target_sources(lahuta_py PRIVATE $<TARGET_OBJECTS:lahuta>)
    get_target_property(_lahuta_libs lahuta LINK_LIBRARIES)
    set_property(TARGET lahuta_py APPEND PROPERTY LINK_LIBRARIES ${_lahuta_libs})
  else()
    target_link_libraries(lahuta_py PRIVATE lahuta)
  endif()
  set_property(TARGET lahuta_py PROPERTY OUTPUT_NAME lahuta)
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



# Set include directories
# include_directories(
#     ${CMAKE_CURRENT_SOURCE_DIR}/include
#     ${CMAKE_CURRENT_SOURCE_DIR}/bond_table
#     ${CMAKE_CURRENT_SOURCE_DIR}/ob
#     ${CMAKE_CURRENT_SOURCE_DIR}/external/ob
#     ${CMAKE_CURRENT_SOURCE_DIR}/external/gemmi/include
#     ${CMAKE_CURRENT_SOURCE_DIR}/external/gemmi/third_party
#     ${CMAKE_CURRENT_SOURCE_DIR}/external/gemmi/third_party/zlib
#     ${RDKit_INCLUDE_DIRS}
#     ${RDKit_ROOT}/include
#     ${EIGEN3_INCLUDE_DIR}
# )
#
# # Collect source files
# file(GLOB_RECURSE SOURCES
#     "src/*.cpp"
#     "src/ob/*.cpp"
#     "src/bond_table/*.cpp"
#     "external/ob/*.cpp"
#     "external/gemmi/src/*.cpp"
# )
# file(GLOB_RECURSE C_SOURCES
#     "external/gemmi/third_party/zlib/*.c"
# )
#
# # Create the executable
# add_executable(${PROJECT_NAME} main.cpp ${SOURCES} ${C_SOURCES})
#
# # Check the option and add the definition for the target
# if (DYNAMIC_CRC_TABLE)
#     target_compile_definitions(${PROJECT_NAME} PRIVATE DYNAMIC_CRC_TABLE)
#     message(STATUS "DYNAMIC_CRC_TABLE is defined")
# else()
#     message(STATUS "DYNAMIC_CRC_TABLE is not defined")
# endif()
#
#
# # Set compile options
# target_compile_options(${PROJECT_NAME} PRIVATE -O3)
#
# Link libraries
# target_link_libraries(${PROJECT_NAME} PRIVATE
# target_link_libraries(lahuta_cpp PRIVATE
#     $ENV{CONDA_PREFIX}/lib/libRDKitGraphMol.so
#     $ENV{CONDA_PREFIX}/lib/libRDKitRDGeneral.so
#     $ENV{CONDA_PREFIX}/lib/libRDKitRDGeometryLib.so
#     $ENV{CONDA_PREFIX}/lib/libRDKitSubstructMatch.so
#     $ENV{CONDA_PREFIX}/lib/libRDKitSmilesParse.so
# )

# Set RPATH for shared libraries
# set_target_properties(lahuta_cpp PROPERTIES
#     INSTALL_RPATH "${RDKit_ROOT}/lib"
#     BUILD_WITH_INSTALL_RPATH TRUE
# )
#
# # Install target
# install(TARGETS lahuta_cpp
#     RUNTIME DESTINATION bin
# )
