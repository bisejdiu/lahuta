function(fetch_gmx_test_data)
    set(SIMDB_DEST "${CMAKE_CURRENT_SOURCE_DIR}/../data/simulationdatabase")
    message(STATUS "Destination for test data: ${SIMDB_DEST}")

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

    foreach(i RANGE ${last_index})
        list(GET GMX_FILENAMES ${i} filename)
        list(GET GMX_URLS ${i} file_url)

        set(dest_file "${SIMDB_DEST}/${filename}")

        if(NOT EXISTS "${dest_file}")
            message(STATUS "Downloading ${filename}...")
            file(DOWNLOAD "${file_url}" "${dest_file}" SHOW_PROGRESS)
        else()
            message(STATUS "Skipping download of ${filename}, file already exists.")
        endif()
    endforeach()
endfunction()
