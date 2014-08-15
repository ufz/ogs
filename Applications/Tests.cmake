
add_test(NAME ogs_no_args COMMAND ogs)
set_tests_properties(ogs_no_args PROPERTIES WILL_FAIL TRUE)

AddOgsBenchmark(ogs_empty_project EmptyProject.xml)
