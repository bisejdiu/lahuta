cmake_minimum_required(VERSION 3.15...3.22)
project(distopia)

set(HWY_ENABLE_TESTS OFF CACHE BOOL "Enable HWY tests" FORCE)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

add_library(libdistopia STATIC "src/distopia.cpp")
cmake_policy(SET CMP0076 NEW)
# target_sources(libdistopia PRIVATE "src/distopia.cpp")
# needed for M_PI symbol on Win
target_compile_definitions(libdistopia PUBLIC -D_USE_MATH_DEFINES)
set_property(TARGET libdistopia PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(libdistopia PUBLIC
        "highway"
        ${CMAKE_CURRENT_LIST_DIR}
        "${CMAKE_CURRENT_LIST_DIR}/include")

message(STATUS "Using vendored HWY subpackage")
# Force Highway to build as static library
set(HWY_FORCE_STATIC_LIBS ON CACHE BOOL "Force Highway to build static" FORCE)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_subdirectory("highway")
# Linking highway publicly so it propagates to the final executable when linking statically
target_link_libraries(libdistopia PUBLIC hwy)

# cmake will append the name lib already , avoid liblibdistopia
set_target_properties(libdistopia PROPERTIES OUTPUT_NAME distopia)
