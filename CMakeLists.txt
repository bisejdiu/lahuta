cmake_minimum_required(VERSION 3.15...3.22)
project(lahuta VERSION 1.0.0 LANGUAGES CXX C)

# Set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add compile definitions
add_compile_definitions(DYNAMIC_CRC_TABLE)

# Set default RDKit_ROOT if not provided
if(NOT DEFINED RDKit_ROOT)
    if(DEFINED ENV{CONDA_PREFIX})
        set(RDKit_ROOT $ENV{CONDA_PREFIX})
    else()
        message(FATAL_ERROR "RDKit_ROOT not set and CONDA_PREFIX environment variable not found. Please set RDKit_ROOT manually.")
    endif()
endif()

# Add RDKit_ROOT to CMAKE_PREFIX_PATH
list(APPEND CMAKE_PREFIX_PATH ${RDKit_ROOT})

# Find external dependencies
find_package(Eigen3 REQUIRED)
find_package(RDKit REQUIRED)

# Set include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/external/gemmi/include
    ${CMAKE_CURRENT_SOURCE_DIR}/bond_table
    ${CMAKE_CURRENT_SOURCE_DIR}/external/gemmi/third_party
    ${CMAKE_CURRENT_SOURCE_DIR}/external/gemmi/third_party/zlib
    ${RDKit_INCLUDE_DIRS}
    ${RDKit_ROOT}/include
    ${EIGEN3_INCLUDE_DIR}
)

# Collect source files
file(GLOB_RECURSE SOURCES
    "src/*.cpp"
    "src/bond_table/*.cpp"
    "external/gemmi/src/*.cpp"
    "external/gemmi/third_party/zlib/*.c"
)

# Create the executable
add_executable(${PROJECT_NAME} main.cpp ${SOURCES})

# Set compile options
target_compile_options(${PROJECT_NAME} PRIVATE -O3)

# Link libraries
target_link_libraries(${PROJECT_NAME} PRIVATE
    $ENV{CONDA_PREFIX}/lib/libRDKitGraphMol.so
    $ENV{CONDA_PREFIX}/lib/libRDKitRDGeneral.so
    $ENV{CONDA_PREFIX}/lib/libRDKitRDGeometryLib.so
    $ENV{CONDA_PREFIX}/lib/libRDKitSubstructMatch.so
    $ENV{CONDA_PREFIX}/lib/libRDKitSmilesParse.so
)

# Set RPATH for shared libraries
set_target_properties(${PROJECT_NAME} PROPERTIES
    INSTALL_RPATH "${RDKit_ROOT}/lib"
    BUILD_WITH_INSTALL_RPATH TRUE
)

# Install target
install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION bin
)
