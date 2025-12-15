#!/bin/bash

set -eu
set -o pipefail
shopt -s inherit_errexit

# where the testee is located
scripts_dev_perf_dir="$(realpath "`dirname "$0"`/..")"

testee="$scripts_dev_perf_dir/util.sh"

ogs_perf_branches=""

. "$scripts_dev_perf_dir/test/helpers.sh"
. "$testee"

test_ogs_perf_cpu_list()
{
    res="`RUN ogs_perf_cpu_list 1`"
    CHECK [ 0 = "$res" ]

    res="`RUN ogs_perf_cpu_list 2`"
    CHECK [ "0,1" = "$res" ]

    res="`RUN ogs_perf_cpu_list 5`"
    CHECK [ "0,1,2,3,4" = "$res" ]
}

test_ogs_perf_conditional_sleep()
{
    res="`RUN ogs_perf_conditional_sleep 1 2 4`"
    CHECK [ -z "$res" ]  # duration < max_no_sleep

    res="`RUN ogs_perf_conditional_sleep 2 2 4`"
    CHECK [ -z "$res" ]  # duration == max_no_sleep

    res="`RUN ogs_perf_conditional_sleep 3 2 8`"
    CHECK [ "waiting for 3 s before continuing" = "$res" ]  # duration > max_no_sleep

    res="`RUN ogs_perf_conditional_sleep 40 2 3`"
    CHECK [ "waiting for 3 s before continuing" = "$res" ]  # duration >> max_no_sleep
}

test_ogs_perf_normalize_revision()
{
    {
        while read input expected_output
        do
            actual_output="`RUN ogs_perf_normalize_revision "$input"`"
            CHECK [ "x$actual_output" = "x$expected_output" ]
        done
    } <<EOF
HEAD HEAD
8417108b738a9da3 8417108b738a9da3
main main
my-remote/my_branch.5 my-remote_my_branch.5
main^ main_
main@{u} main_u_
main~5 main_5
main@{2} main_2_
EOF
}

TEST test_ogs_perf_cpu_list
TEST test_ogs_perf_conditional_sleep
TEST test_ogs_perf_normalize_revision
