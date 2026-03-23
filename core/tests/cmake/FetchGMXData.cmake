function(fetch_gmx_test_data)
    set(SIMDB_DEST "${CMAKE_CURRENT_SOURCE_DIR}/../data/simulationdatabase")

    if(NOT EXISTS "${SIMDB_DEST}")
        file(MAKE_DIRECTORY "${SIMDB_DEST}")
    endif()

    set(GMX_RAW_BASE_URL "https://raw.githubusercontent.com/gromacs/gromacs/8701f03cee72edc28e46b959dacadb2e5b9cc1b4/src/testutils/simulationdatabase")
    set(GMX_BLOB_BASE_URL "https://github.com/gromacs/gromacs/raw/8701f03cee72edc28e46b959dacadb2e5b9cc1b4/src/testutils/simulationdatabase")

    set(GMX_FILENAMES
        "lysozyme.gro"
        "lysozyme.xtc"
        "msd_coords.gro"
        "msd_traj.xtc"
        "msd_traj_rounding_fail.xtc"
        "spc2-traj.gro"
        "spc2-traj.xtc"
    )
    set(GMX_URLS
        "${GMX_RAW_BASE_URL}/lysozyme.gro"
        "${GMX_BLOB_BASE_URL}/lysozyme.xtc"
        "${GMX_RAW_BASE_URL}/msd_coords.gro"
        "${GMX_BLOB_BASE_URL}/msd_traj.xtc"
        "${GMX_BLOB_BASE_URL}/msd_traj_rounding_fail.xtc"
        "${GMX_RAW_BASE_URL}/spc2-traj.gro"
        "${GMX_BLOB_BASE_URL}/spc2-traj.xtc"
    )

    list(LENGTH GMX_FILENAMES num_files)
    math(EXPR last_index "${num_files} - 1")
    set(_lahuta_missing_indices)

    foreach(i RANGE ${last_index})
        list(GET GMX_FILENAMES ${i} filename)
        set(dest_file "${SIMDB_DEST}/${filename}")

        if(NOT EXISTS "${dest_file}")
            list(APPEND _lahuta_missing_indices ${i})
        endif()
    endforeach()

    list(LENGTH _lahuta_missing_indices _lahuta_missing_count)
    if(_lahuta_missing_count EQUAL 0)
        message(STATUS "Using cached GROMACS test data at ${SIMDB_DEST}")
        return()
    endif()

    message(STATUS "Fetching ${_lahuta_missing_count} GROMACS test data file(s) into ${SIMDB_DEST}")

    foreach(i IN LISTS _lahuta_missing_indices)
        list(GET GMX_FILENAMES ${i} filename)
        list(GET GMX_URLS ${i} file_url)
        file(DOWNLOAD "${file_url}" "${SIMDB_DEST}/${filename}")
    endforeach()
endfunction()
