cmake_minimum_required(VERSION 3.20)

project(lahuta LANGUAGES C CXX)
include(CTest)

option(LAHUTA_BUILD_CLI            "Build Lahuta CLI tools" ON)
option(LAHUTA_BUILD_SHARED_CORE    "Build lahuta_core as a shared library" ON)
option(LAHUTA_BUILD_PYTHON         "Build Python bindings" ON)
option(LAHUTA_BUILD_EXAMPLES       "Build Lahuta example executables" OFF)
option(LAHUTA_GENERATE_PY_STUBS    "Generate Python type stubs with pybind11" ON)
option(LAHUTA_VENDOR_BOOST         "Use vendored Boost headers instead of system Boost" ON)
option(LAHUTA_VENDOR_EIGEN         "Use vendored Eigen headers instead of system Eigen" ON)
option(LAHUTA_ENABLE_EXPENSIVE_TESTS "Enable expensive tests that require large datasets" OFF)
option(LAHUTA_ENABLE_LTO           "Enable interprocedural optimization (LTO/IPO) for non-Debug builds" ON)
option(LAHUTA_NATIVE_ARCH          "Enable -march=native for non-Debug builds" OFF)
option(LAHUTA_DISTOPIA_ENABLE_SVE  "Enable Highway SVE targets for distopia" OFF)

add_subdirectory(core)

if(LAHUTA_BUILD_CLI)
  add_subdirectory(cli)
endif()
