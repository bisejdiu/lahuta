set(EIGEN3_VERSION "3.4.1")
set(EIGEN3_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/eigen)

if(NOT EXISTS "${EIGEN3_INCLUDE_DIR}/Eigen/Core")
  message(FATAL_ERROR "Eigen submodule not found at ${EIGEN3_INCLUDE_DIR}. Initialize git submodules and ensure tag ${EIGEN3_VERSION} is checked out.")
endif()

add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
set_target_properties(Eigen3::Eigen PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES ${EIGEN3_INCLUDE_DIR})
