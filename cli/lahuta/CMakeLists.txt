set(lahuta_sources
  main.cpp
  cli/arg_validation.cpp
  cli/global_options.cpp
  commands/contacts.cpp
  commands/createdb.cpp
  commands/extract.cpp
  commands/positions.cpp
)

add_executable(lahuta ${lahuta_sources})

target_include_directories(lahuta PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(lahuta PRIVATE lahuta_core)
target_compile_features(lahuta PRIVATE cxx_std_17)

if(BUILD_TESTING)
  add_test(NAME lahuta_help COMMAND lahuta --help)
  set_tests_properties(lahuta_help PROPERTIES PASS_REGULAR_EXPRESSION "Usage: lahuta \\[global-options\\]")

  add_test(NAME lahuta_contacts_help COMMAND lahuta contacts --help)
  set_tests_properties(lahuta_contacts_help PROPERTIES PASS_REGULAR_EXPRESSION "Compute inter-atomic contacts")

  add_test(NAME lahuta_createdb_help COMMAND lahuta createdb --help)
  set_tests_properties(lahuta_createdb_help PROPERTIES PASS_REGULAR_EXPRESSION "Create a database from alphafold2 models.")

  add_test(NAME lahuta_extract_help COMMAND lahuta extract --help)
  set_tests_properties(lahuta_extract_help PROPERTIES PASS_REGULAR_EXPRESSION "Usage: lahuta extract")

  add_test(NAME lahuta_extract_sequence_help COMMAND lahuta extract --fields sequence --help)
  set_tests_properties(lahuta_extract_sequence_help PROPERTIES PASS_REGULAR_EXPRESSION "Extract")

  add_test(NAME lahuta_positions_help COMMAND lahuta positions --help)
  set_tests_properties(lahuta_positions_help PROPERTIES PASS_REGULAR_EXPRESSION "Usage: lahuta positions")
endif()
