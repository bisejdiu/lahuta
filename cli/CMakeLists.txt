set(lahuta_sources
  main.cpp
  parsing/arg_validation.cpp
  parsing/parsed_args.cpp
  schemas/option_schema.cpp
  schemas/shared_options.cpp
  runner/app.cpp
  runner/pipeline_plan.cpp
  runner/command_runner.cpp
  specs/contacts_spec.cpp
  specs/compaction_rg_spec.cpp
  specs/shape_metrics_spec.cpp
  specs/createdb_spec.cpp
  specs/dssp_spec.cpp
  specs/extract_spec.cpp
  specs/positions_spec.cpp
  specs/quality_metrics_spec.cpp
  specs/sasa_sr_spec.cpp
  tasks/compaction_rg_task.cpp
  tasks/compaction_rg_summary_sink.cpp
  tasks/shape_metrics_task.cpp
  tasks/shape_metrics_summary_sink.cpp
  tasks/contacts_md_source.cpp
  tasks/positions_task.cpp
  tasks/quality_metrics_task.cpp
)

add_executable(lahuta ${lahuta_sources})

target_include_directories(lahuta PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(lahuta PRIVATE lahuta_core)
target_compile_features(lahuta PRIVATE cxx_std_17)

# Set rpath so CLI can find liblahuta_core at runtime.
# In the build tree the executable lives in build/cli and the shared library in
# build/core. Installed layouts differ between wheel and conda/system packages.
if(APPLE)
  set(_cli_build_rpath "@loader_path/../core")
  if(DEFINED SKBUILD)
    set(_cli_install_rpath "@loader_path/../lahuta/lib")
  else()
    set(_cli_install_rpath "@loader_path/../lib")
  endif()
else()
  set(_cli_build_rpath "$ORIGIN/../core")
  if(DEFINED SKBUILD)
    set(_cli_install_rpath "$ORIGIN/../lahuta/lib")
  else()
    set(_cli_install_rpath "$ORIGIN/../lib")
  endif()
endif()
set_target_properties(lahuta PROPERTIES
  BUILD_RPATH "${_cli_build_rpath}"
  INSTALL_RPATH "${_cli_install_rpath}"
  BUILD_WITH_INSTALL_RPATH FALSE
)

install(TARGETS lahuta
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  COMPONENT cli
)

if(BUILD_TESTING)
  add_test(NAME lahuta_help COMMAND lahuta --help)
  set_tests_properties(lahuta_help PROPERTIES PASS_REGULAR_EXPRESSION "Usage: lahuta \\[global-options\\]")

  add_test(NAME lahuta_contacts_help COMMAND lahuta contacts --help)
  set_tests_properties(lahuta_contacts_help PROPERTIES PASS_REGULAR_EXPRESSION "Compute inter-atomic contacts")

  add_test(NAME lahuta_createdb_help COMMAND lahuta createdb --help)
  set_tests_properties(lahuta_createdb_help PROPERTIES PASS_REGULAR_EXPRESSION "Create a database from alphafold2 models.")

  add_test(NAME lahuta_quality_metrics_help COMMAND lahuta quality-metrics --help)
  set_tests_properties(lahuta_quality_metrics_help PROPERTIES PASS_REGULAR_EXPRESSION "Usage: lahuta quality-metrics")

  add_test(NAME lahuta_extract_help COMMAND lahuta extract --help)
  set_tests_properties(lahuta_extract_help PROPERTIES PASS_REGULAR_EXPRESSION "Usage: lahuta extract")

  add_test(NAME lahuta_extract_sequence_help COMMAND lahuta extract --fields sequence --help)
  set_tests_properties(lahuta_extract_sequence_help PROPERTIES PASS_REGULAR_EXPRESSION "Extract")

  add_test(NAME lahuta_positions_help COMMAND lahuta positions --help)
  set_tests_properties(lahuta_positions_help PROPERTIES PASS_REGULAR_EXPRESSION "Usage: lahuta positions")

  add_test(NAME lahuta_compaction_rg_help COMMAND lahuta compaction-rg --help)
  set_tests_properties(lahuta_compaction_rg_help PROPERTIES PASS_REGULAR_EXPRESSION "Usage: lahuta compaction-rg")

  add_test(NAME lahuta_shape_metrics_help COMMAND lahuta shape-metrics --help)
  set_tests_properties(lahuta_shape_metrics_help PROPERTIES PASS_REGULAR_EXPRESSION "Usage: lahuta shape-metrics")

  add_test(NAME lahuta_sasa_sr_help COMMAND lahuta sasa-sr --help)
  set_tests_properties(lahuta_sasa_sr_help PROPERTIES PASS_REGULAR_EXPRESSION "Usage: lahuta sasa-sr")
endif()
