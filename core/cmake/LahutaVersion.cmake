#
# Lahuta version management
# Reads VERSION file from repository root and propagates to C++ and Python.
#
function(lahuta_read_version_file)
  set(LAHUTA_VERSION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/../VERSION")

  if(NOT EXISTS "${LAHUTA_VERSION_FILE}")
    message(FATAL_ERROR "Lahuta version file not found at ${LAHUTA_VERSION_FILE}.")
  endif()

  file(READ "${LAHUTA_VERSION_FILE}" _lahuta_version_raw)
  string(STRIP "${_lahuta_version_raw}" LAHUTA_VERSION_FULL)

  if(NOT LAHUTA_VERSION_FULL)
    message(FATAL_ERROR "Lahuta version file ${LAHUTA_VERSION_FILE} is empty.")
  endif()

  # Parse semantic version components
  string(REGEX MATCH "^([0-9]+)\\.([0-9]+)\\.([0-9]+)" _ "${LAHUTA_VERSION_FULL}")
  if(NOT CMAKE_MATCH_COUNT)
    message(FATAL_ERROR "Lahuta version '${LAHUTA_VERSION_FULL}' must start with semantic version 'major.minor.patch'.")
  endif()

  # Export to parent scope
  set(LAHUTA_VERSION_MAJOR "${CMAKE_MATCH_1}" PARENT_SCOPE)
  set(LAHUTA_VERSION_MINOR "${CMAKE_MATCH_2}" PARENT_SCOPE)
  set(LAHUTA_VERSION_PATCH "${CMAKE_MATCH_3}" PARENT_SCOPE)
  set(LAHUTA_VERSION_BASE  "${CMAKE_MATCH_0}" PARENT_SCOPE)

  # Extract any suffix (e.g., "-dev", "-alpha", etc.)
  string(REGEX REPLACE "^${CMAKE_MATCH_0}" "" LAHUTA_VERSION_SUFFIX "${LAHUTA_VERSION_FULL}")
  set(LAHUTA_VERSION_SUFFIX "${LAHUTA_VERSION_SUFFIX}" PARENT_SCOPE)
  set(LAHUTA_VERSION_FULL "${LAHUTA_VERSION_FULL}" PARENT_SCOPE)

  message(STATUS "Configuring Lahuta ${LAHUTA_VERSION_FULL}")
endfunction()

# Configure version information for the C++ core library target
function(lahuta_configure_core_version target_name)
  if(NOT TARGET ${target_name})
    message(FATAL_ERROR "Target '${target_name}' does not exist.")
  endif()

  set_target_properties(${target_name} PROPERTIES
    VERSION ${LAHUTA_VERSION_BASE}
    SOVERSION ${LAHUTA_VERSION_MAJOR}
  )

  target_compile_definitions(${target_name} PUBLIC
    LAHUTA_VERSION_MAJOR=${LAHUTA_VERSION_MAJOR}
    LAHUTA_VERSION_MINOR=${LAHUTA_VERSION_MINOR}
    LAHUTA_VERSION_PATCH=${LAHUTA_VERSION_PATCH}
    LAHUTA_VERSION_SUFFIX="${LAHUTA_VERSION_SUFFIX}"
    LAHUTA_VERSION_STRING="${LAHUTA_VERSION_FULL}"
  )
endfunction()

# Configure Python _version.py file for direct CMake installs
function(lahuta_configure_python_version py_pkg_root output_dir)
  set(version_template "${py_pkg_root}/lahuta/_version.py")
  set(version_output "${output_dir}/_version.py")

  if(NOT EXISTS "${version_template}")
    message(FATAL_ERROR "Python version template not found at ${version_template}")
  endif()

  set(LAHUTA_CMAKE_CONFIGURED "TRUE")
  configure_file("${version_template}" "${version_output}" @ONLY)

  # Export the output path to parent scope so it can be installed
  set(LAHUTA_PYTHON_VERSION_FILE "${version_output}" PARENT_SCOPE)
endfunction()
