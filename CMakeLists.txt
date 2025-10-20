cmake_minimum_required(VERSION 3.20)

project(lahuta LANGUAGES C CXX)
include(CTest)

option(LAHUTA_BUILD_CLI "Build Lahuta CLI tools" ON)

add_subdirectory(core)

if(LAHUTA_BUILD_CLI)
  add_subdirectory(cli)
endif()
