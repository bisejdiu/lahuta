cmake_minimum_required(VERSION 3.15...3.22)
project(distopia)

set(HWY_ENABLE_TESTS OFF CACHE BOOL "Enable HWY tests" FORCE)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

add_library(libdistopia SHARED "src/distopia.cpp")
cmake_policy(SET CMP0076 NEW)
# target_sources(libdistopia PRIVATE "src/distopia.cpp")
# needed for M_PI symbol on Win
target_compile_definitions(libdistopia PUBLIC -D_USE_MATH_DEFINES)
target_include_directories(libdistopia PUBLIC
        "highway"
        ${CMAKE_CURRENT_LIST_DIR}
        "${CMAKE_CURRENT_LIST_DIR}/include")

message(STATUS "Using vendored HWY subpackage")
add_subdirectory("highway")
target_link_libraries(libdistopia PRIVATE hwy)

# cmake will append the name lib already , avoid liblibdistopia
set_target_properties(libdistopia PROPERTIES OUTPUT_NAME distopia)
if(APPLE)
    set_target_properties(libdistopia PROPERTIES INSTALL_RPATH "@loader_path")
else()
    set_target_properties(libdistopia PROPERTIES INSTALL_RPATH "\$ORIGIN")
endif()
