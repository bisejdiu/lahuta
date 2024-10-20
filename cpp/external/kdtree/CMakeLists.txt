# external/kdtree/CMakeLists.txt
cmake_minimum_required(VERSION 3.6)
project(KDTree CXX)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

add_library(kdtree KDTree.cpp)

target_include_directories(kdtree PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_options(kdtree PRIVATE
    $<$<CONFIG:RELEASE>:-O3 -DNDEBUG>
)
