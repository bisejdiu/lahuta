set(LAHUTA_EXAMPLE_SOURCES
  getcontacts_example.cpp
  md_pipeline_example.cpp
  md_stage_manager_example.cpp
  nmr_pipeline_example.cpp
  nmr_stage_manager_example.cpp
  use_resolver.cpp
)

foreach(example_src IN LISTS LAHUTA_EXAMPLE_SOURCES)
  get_filename_component(example_target "${example_src}" NAME_WE)
  add_executable("${example_target}" "${example_src}")
  target_link_libraries("${example_target}" PRIVATE lahuta_core)
  target_compile_features("${example_target}" PRIVATE cxx_std_17)
endforeach()
