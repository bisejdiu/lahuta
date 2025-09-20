if(NOT BUILD_TESTING)
  return()
endif()

include(FetchContent)

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

add_executable(lahuta_tests ${LAHUTA_TEST_MAIN} ${LAHUTA_TEST_UNITS})
target_compile_features(lahuta_tests PRIVATE cxx_std_17)

target_include_directories(lahuta_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if(MSVC)
  target_compile_options(lahuta_tests PRIVATE /fp:precise)
else()
  target_compile_options(lahuta_tests PRIVATE -fno-fast-math -fexcess-precision=standard -ffp-contract=off)
endif()

target_link_libraries(lahuta_tests
  PRIVATE
    lahuta_core
    GTest::gtest
    ${_gmock_target}
)

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
