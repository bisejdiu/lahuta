set(LAHUTA_PROG_SOURCES
  get_organism_data.cpp
  get_sequence_data.cpp
  get_plddt_data.cpp
  get_dssp_data.cpp
  get_position_data.cpp
)

foreach(prog_src IN LISTS LAHUTA_PROG_SOURCES)
  get_filename_component(prog_target "${prog_src}" NAME_WE)
  add_executable("${prog_target}" "${prog_src}")
  target_link_libraries("${prog_target}" PRIVATE lahuta_core)
  target_compile_features("${prog_target}" PRIVATE cxx_std_17)
endforeach()
