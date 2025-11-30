#!/bin/sh

set -eu
set -o pipefail
shopt -s inherit_errexit

exists() { [ -e "$1" ]; }

### some definitions
#
today="`date +%F`"

# to make dirs/branches unique
uuid="`uuidgen`"
#
# temporary branch
branch="test-master-$today-$uuid"

# temporary git worktree
worktree="../ogs-ogs-perf/tmpworktree-$uuid"

# where the testee is located
scripts_dev_perf_dir="$(realpath "`dirname "$0"`/..")"

. "$scripts_dev_perf_dir/test/helpers.sh"

testee="$scripts_dev_perf_dir/run_ctest.sh"

# compiler name occurs in generated dir names
comp=gcc


### test definitions
setup()
{
    # consistency check: measurement results dir must not exist
    CHECK_NOT exists "measurements-$today"*/

    # make sure the master branch exists, which it might not in CI jobs
    git fetch origin master

    mkdir -p ../ogs-ogs-perf  # would be generated anyways by testee

    CHECK_NOT [ -e "$worktree" ]
    git worktree add --checkout -b "$branch" "$worktree" origin/master  # prepare a worktree and create a new branch

    # determine commit hash
    hash="`git rev-parse "$branch"`"

    # check that we are in the "ogs" dir
    CHECK [ ogs = "$(basename "`pwd`")" ]

    CHECK [ -n "$CPM_SOURCE_CACHE" ]  # env var must be set

    # CI var presence changes where CTest log files are put
    unset CI
}


# only called if all tests ran successfully
teardown()
{
    # clean up to simplify test development
    # didn't exist before (checked in setup), exists now
    rm -r "measurements-$today"*/

    git worktree remove "$worktree"
    CHECK_NOT [ -e "$worktree" ]

    git branch -D "$branch"
    git branch -D "$branch-2"

    rm -r "../ogs-ogs-perf/$today-$branch-$comp"
}

test_0_help()
{
    msg="`RUN "$testee" --help | head -n1`"
    CHECK [ "Usage: run_ctest.sh [CTEST_ARGS...]" = "$msg" ]
}

test_1_list_only()
{
    RUN env \
        ogs_perf_branches="$branch" ogs_perf_num_test_runs=1 \
        "$testee" -N

    CHECK [ -d ../ogs-ogs-perf ]  # perf dir must exist
    CHECK [ -d "../ogs-ogs-perf/$today-$branch-$comp/ogs" ]
    CHECK [ -d "../ogs-ogs-perf/$today-$branch-$comp/build/release" ]
    CHECK [ -d "measurements-$today" ]  # measurements dir must exist
    CHECK [ -e "measurements-$today/$today-$branch-$hash/logs/ctest-$branch.log" ]  # ctest stdout log file must exist
    CHECK_NOT exists "measurements-$today/$today-$branch-$hash/logs."*  # single run - no other log dir must exist
}

test_2_run_some_tests()
{
    RUN env \
        ogs_perf_branches="$branch" ogs_perf_num_test_runs=3 \
        "$testee" -R ogs-Mechanics/Linear/square

    # ctest stdout files
    CHECK [ -e "measurements-$today/$today-$branch-$hash/logs/ctest-$branch.log" ]
    CHECK [ -e "measurements-$today/$today-$branch-$hash/logs.~1~/ctest-$branch.log" ]
    CHECK [ -e "measurements-$today/$today-$branch-$hash/logs.~2~/ctest-$branch.log" ]
    CHECK [ -e "measurements-$today/$today-$branch-$hash/logs.~3~/ctest-$branch.log" ]

    # OGS log files
    CHECK exists "measurements-$today/$today-$branch-$hash/logs/ogs-Mechanics_Linear_square"*.txt
    CHECK exists "measurements-$today/$today-$branch-$hash/logs.~2~/ogs-Mechanics_Linear_square"*.txt
    CHECK exists "measurements-$today/$today-$branch-$hash/logs.~3~/ogs-Mechanics_Linear_square"*.txt
}

test_3_modify_and_build()
{
    # prerequisites/preparation
    local ccache_statslog="../ogs-ogs-perf/$today-$branch-$comp/build/release_ccache_stats_$today.log"
    CHECK [ -e "$ccache_statslog" ]  # ccache statslog must exist
    rm "$ccache_statslog"  # remove it (zeroes stats)

    # modify a file
    echo "const char test_`echo "$uuid" | tr -d -`[] = \"$uuid\";" >>"$worktree/Applications/CLI/ogs.cpp"
    git -C "$worktree" commit -a -m "test-commit" --no-verify
    git --no-pager -C "$worktree" log -n4 --oneline --color=always

    RUN env \
        ogs_perf_branches="$branch" ogs_perf_num_test_runs=3 \
        "$testee" -R ogs-Mechanics/Linear/square

    # checks
    CHECK [ -e "$ccache_statslog" ]  # ccache statslog must exist again, i.e. ccache must have been called

    # ogs.cpp was modified. There should be 0/1 ccache hits (reported two/three times: calls/local storage/maybe remote storage)
    statslog="$(CCACHE_STATSLOG="`realpath "$ccache_statslog"`" ccache --show-log-stats)"
    echo "$statslog"
    CHECK [ 2 -le "$(echo "$statslog" | grep "Hits: *0 / 1" | wc -l)" ]
    CHECK [ 2 -le "$(echo "$statslog" | grep "Misses: *1 / 1" | wc -l)" ]
}

test_4_revert_and_build()
{
    # prerequisites/preparation
    local ccache_statslog="../ogs-ogs-perf/$today-$branch-$comp/build/release_ccache_stats_$today.log"
    rm "$ccache_statslog"  # remove statslog (zeroes stats)

    git -C "$worktree" revert HEAD --no-edit
    git --no-pager -C "$worktree" log -n4 --oneline --color=always

    # call testee
    RUN env \
        ogs_perf_branches="$branch" ogs_perf_num_test_runs=3 \
        "$testee" -R ogs-Mechanics/Linear/square

    # checks
    CHECK [ -e "$ccache_statslog" ]  # ccache statslog must exist again, i.e. ccache must have been called
    # ogs.cpp was reverted. There should be 1/1 ccache hits (reported two times)
    local statslog="$(CCACHE_STATSLOG="`realpath "$ccache_statslog"`" ccache --show-log-stats)"
    echo "$statslog"
    CHECK [ 2 -eq "$(echo "$statslog" | grep "Hits: *1 / 1" | wc -l)" ]
    CHECK [ 2 -eq "$(echo "$statslog" | grep "Misses: *0 / 1" | wc -l)" ]
}

test_5_build_different_branch_name()
{
    # prerequisite/preparation
    git -C "$worktree" branch -c "$branch" "$branch-2"

    RUN env \
        ogs_perf_branches="$branch-2" ogs_perf_num_test_runs=1 \
        "$testee" -R ogs-Mechanics/Linear/square

    # checks
    ccache_statslog2="../ogs-ogs-perf/$today-$branch-2-$comp/build/release_ccache_stats_$today.log"
    CHECK [ -e "$ccache_statslog2" ]  # ccache statslog must exist
    statslog="$(CCACHE_STATSLOG="`realpath "$ccache_statslog2"`" ccache --show-log-stats)"
    echo "$statslog"
    # expect more than 90% ccache hits
    CHECK [ 90 -lt "$(echo "$statslog" | sed -n -e "/Hits:/{s/^ *Hits:.*( *\([0-9]*\)\..*%)/\1/p; q}")" ]

    # won't be reused, remove
    rm -r "../ogs-ogs-perf/$today-$branch-2-$comp"
}

test_6_remote_tracking_branch()
{
    local b="origin/master"
    local B="origin_master" # normalized name
    RUN env \
        ogs_perf_branches="$b" ogs_perf_num_test_runs=1 \
        "$testee" -N -R ogs-Mechanics

    CHECK [ -d ../ogs-ogs-perf ]  # perf dir must exist
    CHECK [ -d "../ogs-ogs-perf/$today-$B-$comp/ogs" ]
    CHECK [ -d "../ogs-ogs-perf/$today-$B-$comp/build/release" ]
    CHECK [ -d "measurements-$today" ]  # measurements dir must exist
    CHECK [ -e "measurements-$today/$today-$B-$hash/logs/ctest-$B.log" ]  # ctest stdout log file must exist

    # fixed branch name in $b, multiple log dirs might exist (e.g., from
    # previous test runs), hence test disabled
    # CHECK_NOT exists "measurements-$today/$today-$b-$hash/logs."*

    # won't be reused, remove
    rm -r "../ogs-ogs-perf/$today-$B-$comp"
}

test_7_caret_ref()
{
    # should be run after modifications in test 3 or 4 to guarantee that the
    # selected commit can be built
    local b="$branch^"
    local B="${branch}_" # normalized name
    local h="`git rev-parse "$b"`"
    RUN env \
        ogs_perf_branches="$b" ogs_perf_num_test_runs=1 \
        "$testee" -N -R ogs-Mechanics

    CHECK [ -d ../ogs-ogs-perf ]  # perf dir must exist
    CHECK [ -d "../ogs-ogs-perf/$today-$B-$comp/ogs" ]
    CHECK [ -d "../ogs-ogs-perf/$today-$B-$comp/build/release" ]
    CHECK [ -d "measurements-$today" ]  # measurements dir must exist
    CHECK [ -e "measurements-$today/$today-$B-$h/logs/ctest-$B.log" ]  # ctest stdout log file must exist
    CHECK_NOT exists "measurements-$today/$today-$B-$h/logs."*  # single run - no other log dir must exist

    # won't be reused, remove
    rm -r "../ogs-ogs-perf/$today-$B-$comp"
}

test_8_commit_hash()
{
    # should be run after modifications in test 3 or 4 to guarantee that the
    # selected commit hash is unique for this test run (and not the same as
    # origin/master)
    local h="`git rev-parse "$branch"`"
    local b="$h"
    RUN env \
        ogs_perf_branches="$b" ogs_perf_num_test_runs=1 \
        "$testee" -N -R ogs-Mechanics

    CHECK [ -d ../ogs-ogs-perf ]  # perf dir must exist
    CHECK [ -d "../ogs-ogs-perf/$today-$b-$comp/ogs" ]
    CHECK [ -d "../ogs-ogs-perf/$today-$b-$comp/build/release" ]
    CHECK [ -d "measurements-$today" ]  # measurements dir must exist
    CHECK [ -e "measurements-$today/$today-$b-$h/logs/ctest-$b.log" ]  # ctest stdout log file must exist
    CHECK_NOT exists "measurements-$today/$today-$b-$h/logs."*  # single run - no other log dir must exist

    # won't be reused, remove
    rm -r "../ogs-ogs-perf/$today-$b-$comp"
}


RUN setup

TEST test_0_help

TEST test_1_list_only
TEST test_2_run_some_tests

# the following tests check that ccache is effective
TEST test_3_modify_and_build
TEST test_4_revert_and_build
TEST test_5_build_different_branch_name

TEST test_6_remote_tracking_branch
TEST test_7_caret_ref
TEST test_8_commit_hash

RUN teardown
