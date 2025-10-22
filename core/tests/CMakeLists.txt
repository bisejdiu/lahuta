if(NOT BUILD_TESTING)
  return()
endif()

include(FetchContent)

option(FETCH_TESTDATA "Fetch external testdata (GROMACS simulationdatabase) when building tests" ON)
if(FETCH_TESTDATA)
  include(cmake/FetchGMXData.cmake)
  fetch_gmx_test_data()
endif()

set(INSTALL_GTEST OFF CACHE BOOL "Disable install for gtest" FORCE)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.17.0
)

FetchContent_MakeAvailable(googletest)

set(LAHUTA_TEST_MAIN main.cpp)
# Collect all test sources recursively. Exclude the main compilation unit to avoid duplicate mains
file(GLOB_RECURSE LAHUTA_TEST_UNITS CONFIGURE_DEPENDS
  *.cpp
)
list(REMOVE_ITEM LAHUTA_TEST_UNITS ${CMAKE_CURRENT_SOURCE_DIR}/${LAHUTA_TEST_MAIN})

# Requires large manual data downloads
if(NOT FETCH_TESTDATA)
  list(REMOVE_ITEM LAHUTA_TEST_UNITS ${CMAKE_CURRENT_SOURCE_DIR}/md/test_gmx_compat_xtc_gro.cpp)
  message(WARNING "FETCH_TESTDATA is disabled. GROMACS compatibility tests will be skipped. Set -DFETCH_TESTDATA=ON to enable them.")
endif()

# Exclude expensive & local-data dependent tests
if(NOT (DEFINED ENV{LAHUTA_ENABLE_EXPENSIVE_TESTS} AND NOT "$ENV{LAHUTA_ENABLE_EXPENSIVE_TESTS}" STREQUAL ""))
  list(REMOVE_ITEM LAHUTA_TEST_UNITS ${CMAKE_CURRENT_SOURCE_DIR}/pipeline/test_topology_cache_validation.cpp)
  list(REMOVE_ITEM LAHUTA_TEST_UNITS ${CMAKE_CURRENT_SOURCE_DIR}/pipeline/test_xtc_parallel_coordinates.cpp)
endif()

add_executable(lahuta_tests ${LAHUTA_TEST_MAIN} ${LAHUTA_TEST_UNITS})
target_compile_features(lahuta_tests PRIVATE cxx_std_17)

target_include_directories(lahuta_tests PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/cli/lahuta
)

if(MSVC)
  target_compile_options(lahuta_tests PRIVATE /fp:precise)
else()
  target_compile_options(lahuta_tests PRIVATE -fno-fast-math -ffp-contract=off)
endif()

target_link_libraries(lahuta_tests
  PRIVATE
    lahuta_core
    GTest::gtest
    ${_gmock_target}
)

if(DEFINED ENV{LAHUTA_ENABLE_EXPENSIVE_TESTS} AND NOT "$ENV{LAHUTA_ENABLE_EXPENSIVE_TESTS}" STREQUAL "")
  set(_lahuta_simdb_dir "${CMAKE_SOURCE_DIR}/gpcrmd")
else()
  set(_lahuta_simdb_dir "${CMAKE_CURRENT_SOURCE_DIR}/data/simulationdatabase")
endif()

target_compile_definitions(lahuta_tests PRIVATE LAHUTA_SIMDB_DIR=\"${_lahuta_simdb_dir}\")

# Sanitizers and coverage
option(ENABLE_ASAN "Enable AddressSanitizer + UndefinedBehaviorSanitizer" OFF)
option(ENABLE_TSAN "Enable ThreadSanitizer" OFF)
option(ENABLE_COVERAGE "Enable Clang source-based coverage" OFF)

# Instrument both the tests and the code-under-test (lahuta_core)
if(ENABLE_ASAN AND NOT ENABLE_TSAN)
  target_compile_options(lahuta_tests PRIVATE -O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined)
  target_link_options(   lahuta_tests PRIVATE -fno-omit-frame-pointer -fsanitize=address,undefined)

  target_compile_options(lahuta_core  PRIVATE -O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined)
  target_link_options(   lahuta_core  PRIVATE -fno-omit-frame-pointer -fsanitize=address,undefined)
endif()

if(ENABLE_TSAN AND NOT ENABLE_ASAN)
  target_compile_options(lahuta_tests PRIVATE -O1 -g -fno-omit-frame-pointer -fsanitize=thread)
  target_link_options(   lahuta_tests PRIVATE -fno-omit-frame-pointer -fsanitize=thread)

  target_compile_options(lahuta_core  PRIVATE -O1 -g -fno-omit-frame-pointer -fsanitize=thread)
  target_link_options(   lahuta_core  PRIVATE -fno-omit-frame-pointer -fsanitize=thread)
endif()

if(ENABLE_COVERAGE)
  # Apply clang source-based coverage to both tests and library
  target_compile_options(lahuta_tests PRIVATE -fprofile-instr-generate -fcoverage-mapping)
  target_link_options(   lahuta_tests PRIVATE -fprofile-instr-generate -fcoverage-mapping)

  target_compile_options(lahuta_core  PRIVATE -fprofile-instr-generate -fcoverage-mapping)
  target_link_options(   lahuta_core  PRIVATE -fprofile-instr-generate -fcoverage-mapping)
endif()

include(GoogleTest)
gtest_discover_tests(lahuta_tests
  DISCOVERY_MODE PRE_TEST
  PROPERTIES TIMEOUT 30
)
