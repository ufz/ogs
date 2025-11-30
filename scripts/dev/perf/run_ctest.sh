#!/bin/sh

script_dir="$(dirname "$(realpath "$0")")"

help()
{
    cat <<EOF
Usage: `basename "$0"` [CTEST_ARGS...]

Clones the OGS source code to a new git repo, builds OGS and runs CTest a number
of times. All commandline arguments are passed on to CTest.

The working directory of this script is expected to be inside an OGS git
repository (typical directory name: ogs). The source and build directories will
be created next to that repository in a subdirectory of an ogs-ogs-perf directory.

This script copies the log files of CTests to a subdirectory of
measurements-<yyyy>-<mm>-<dd> in the current directory.

The behaviour of this script is controlled by some environment variables, see
`dirname "$0"`/settings.sh, most notably:

- ogs_perf_num_test_runs - the number of ctest runs
- ogs_perf_branches - the branches/commits to test (space separated list)
- ogs_perf_compiler - the compiler being used (gcc/clang)
EOF
}

if [ "x$1" = "x--help" ]
then
    help
    exit
fi

. "$script_dir/util.sh"

# use this for debugging
#set -x
#export PS4='+<(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}()> }'

do_ctest()
{
    (
        set -e
        shopt -s inherit_errexit

        for branch in $ogs_perf_branches
        do
            ogs_perf_detect_git_and_build_and_command "$branch"

            echo "---- ctests for branch $branch (hash=$ogs_perf_commit_hash) --------"

            local start="`date +%s`"
            ogs_perf_setup "$ogs_perf_commit_hash" "$ogs_perf_git_root" "$ogs_perf_build_dir"
            local end="`date +%s`"

            local duration="$((end - start))"
            echo "setup took $duration s"

            sync
            ogs_perf_conditional_sleep "$duration" 8 20

            git --no-pager -C "$ogs_perf_git_root" log -n 4 --oneline

            for i in `seq "$ogs_perf_num_test_runs"`;
            do
                echo "---- ctest run #$i/$ogs_perf_num_test_runs for branch $branch --------"
                (
                    set +e
                    cd "$ogs_perf_build_dir"

                    rm -rf logs/
                    mkdir logs

                    local log_file="logs/ctest-$ogs_perf_branch_normalized.log"

                    OMP_NUM_THREADS=1 \
                        taskset -c 0-7 \
                        ctest -E LARGE "$@" 2>&1 | tee "$log_file"
                )

                local log_dir="measurements-$ogs_perf_start_date/$ogs_perf_start_date-$ogs_perf_branch_normalized-$ogs_perf_commit_hash"
                mkdir -p "$log_dir"
                mv "$ogs_perf_build_dir/logs" "$log_dir/" --backup=numbered
                echo "ctest logs have been written to `realpath "$log_dir"`"
            done
        done
    )
}

do_ctest "$@"
