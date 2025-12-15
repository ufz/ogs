# Number of test repetitions.
export ogs_perf_num_test_runs="${ogs_perf_num_test_runs:-3}"

# Space separated list of branches/commit hashes to test.
export ogs_perf_branches="${ogs_perf_branches?-Error: the variable ogs_perf_branches is not set.}"

# Compiler used.
export ogs_perf_compiler="${ogs_perf_compiler:-gcc}"

# Start date of the current test run.
# Needed such that all outputs from a long running test will be put in the same
# directories
export ogs_perf_start_date="${ogs_perf_start_date:-$(date +%F)}"
