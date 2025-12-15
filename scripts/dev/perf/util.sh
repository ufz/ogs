# helpers for ogs performance measurements

script_dir="$(dirname "$(realpath "$BASH_SOURCE")")"
. "$script_dir/settings.sh"


error()
{
    echo "Error: $*" >&2
}

warn()
{
    echo "Warning: $*" >&2
}

die()
{
    error "$*"
    exit 1
}

# Generates a list of comma separated numbers from 0 to n-1.
#
# Arguments:
# - n - length of the generated list
#
# Returns:
# - prints a comma separated list of integers on stdout
#
# The list can be used to be used with taskset -c, hence the name of the
# function.
#
ogs_perf_cpu_list()
{
    local min_cpu=0
    local threads="$1"
    seq -s, $min_cpu $((min_cpu + threads - 1))
}

# Runs a command on a fixed number of CPUs
#
# Arguments:
# - n - how many CPUs to pin the command to (0 ... n-1)
# - command... - the remaining arguments define the command to be run
#
ogs_perf_run_with_fixed_cpus()
{
    local threads="$1"
    shift
    taskset -c "`ogs_perf_cpu_list "$threads"`" "$@"
}

# Pins a running command to a fixed number of CPUs
#
# Arguments:
# - pid - the pid of the command
# - n - how many CPUs to pin the command to (0 ... n-1)
#
ogs_perf_fix_cpus()
{
    local pid="$1"
    local threads="$2"
    taskset -cp "`ogs_perf_cpu_list "$threads"`" "$pid"
}

# Initializes an OGS build directory
#
# Arguments:
# - srcdir - the source directory, must exist
# - builddir - the build directory, must not exist, yet
#
# Special care is taken to set CCache and CPM options that allow for a high
# ratio of CCache hits. Furthermore, compiler options are set to yield both well
# optimized code and enough info for performance assessments (no omission of
# frame pointers, debug info)
#
ogs_perf_init_build_dir()
{
    local srcdir="$1"
    local builddir="$2"
    local preset=release
    local cc=gcc
    local cxx=g++
    if [ "x$ogs_perf_compiler" = "xclang" ]
    then
        cc=clang
        cxx=clang++
    elif [ "x$ogs_perf_compiler" != "xgcc" ]
    then
        die "unsupported \$ogs_perf_compiler: $ogs_perf_compiler"
    fi

    if [ ! -d "$srcdir" ]; then echo "error srcdir '$srcdir' does not exist" >&2; exit 1; fi
    if [ -d "$builddir" ]; then echo "error builddir '$builddir' already exists" >&2; exit 1; fi
    CXXFLAGS="-Wp,-D_GLIBCXX_ASSERTIONS -fno-omit-frame-pointer -mno-omit-leaf-frame-pointer"  \
        cmake -S "$srcdir" -B "$builddir" --preset "$preset" \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        "-DCCACHE_ENV=CCACHE_BASEDIR=`realpath "$srcdir/.."`;CCACHE_SLOPPINESS=pch_defines,time_macros;CCACHE_NOHASHDIR='';CCACHE_STATSLOG=${builddir%/}_ccache_stats_`date +%F`.log" \
        -DBUILD_SHARED_LIBS=ON \
        -DOGS_EIGEN_DYNAMIC_SHAPE_MATRICES=OFF \
        -DOGS_EIGEN_INITIALIZE_MATRICES_BY_NAN=OFF \
        -DOGS_USE_MKL=OFF \
        -DOGS_USE_MFRONT=OFF \
        -DOGS_USE_POETRY=OFF \
        -DOGS_USE_UNITY_BUILDS=OFF \
        -DCMAKE_DISABLE_PRECOMPILE_HEADERS=ON \
        -DCPM_SOURCE_CACHE="$CPM_SOURCE_CACHE" \
        -DCMAKE_C_COMPILER="$cc" -DCMAKE_CXX_COMPILER="$cxx" \
        -DCMAKE_CXX_FLAGS_RELWITHDEBINFO="-O3 -g -DNDEBUG" \
        -DCMAKE_C_FLAGS_RELWITHDEBINFO="-O3 -g -DNDEBUG"
}

# Detects source and build directories and sets certain ogs_perf_... variables
# used by other functions.
#
# Arguments:
# - branch - the branch name or commit hash to be built later on
#
# Sets the following variables:
# - ogs_perf_git_root - the source dir, created as a shared git clone if it does
#                       not exist, yet
# - ogs_perf_build_dir - the build dir to be used
# - ogs_perf_commit_hash - the commit hash of the passed branch
# - ogs_perf_branch_normalized - the branch name with some characters replaced
#                                by underscores; to be used, e.g., in file paths
#
ogs_perf_detect_git_and_build_and_command()
{
    local branch="$1"
    local git_root=""

    # assuming we are in the OGS git repo
    git_root="`git rev-parse --show-toplevel`"
    local base_git="`basename "$git_root"`"

    # note: commit hash is resolved in the current repo (i.e., the OGS git repo)
    # "out arg"
    ogs_perf_commit_hash="`git rev-parse "$branch"`"

    # "out arg"
    ogs_perf_branch_normalized="`ogs_perf_normalize_revision "$branch"`"

    local src_and_build_root="$git_root-ogs-perf/`date +%F`-$ogs_perf_branch_normalized-$ogs_perf_compiler"

    # using directory names ($base_git, build/release) that should be in use
    # during the usual development work anyways -> higher likelihood of
    # cache hits
    #
    # "out arg"
    ogs_perf_git_root="$src_and_build_root/$base_git"
    # "out arg"
    ogs_perf_build_dir="$src_and_build_root/build/release"

    if ! [ -e "$ogs_perf_git_root" ]
    then
        git clone --shared --no-checkout "$git_root" "$ogs_perf_git_root"
        # --detach in case $branch is not a local branch, but, e.g. a remote
        # tracking branch or a commit hash
        git -C "$ogs_perf_git_root" switch "$ogs_perf_commit_hash" -q --detach
    else
        local alt="`cat "$ogs_perf_git_root/.git/objects/info/alternates"`"
        if [ "x$alt" != "x$git_root/.git/objects" ]; then
            # git common dir is used with git worktree, for instance
            local git_com="`realpath "$(git rev-parse --git-common-dir)"`"

            # We only emit a warning message, because, e.g., git does not
            # create a shared clone of a shallow repository. This is the
            # case in the CI jobs.
            [ "x$alt" = "x$git_com/objects" ] \
                || warn "'$ogs_perf_git_root' is not a shared git clone of '$git_root'"
            # note: if a warning is printed, $git_root itself might be a
            # shared clone
        fi

        local url="`git -C "$ogs_perf_git_root" remote get-url origin`"
        [ "x$url" = "x$git_root" ] \
            || die "remote 'origin' does not point to '$git_root' but to '$url'"

        git -C "$ogs_perf_git_root" fetch origin -q
        git -C "$ogs_perf_git_root" switch "$ogs_perf_commit_hash" -q --detach
    fi
}

# Builds OGS, creating the build directory if necessary.
#
# Arguments:
# - commit - the commit (hash) to be built
# - git_root - the source directory
# - build_dir - the build directory
#
# The build log is saved in a file next to the build directory
# (...-ogs-perf-build.log)
#
ogs_perf_setup()
{
    local commit="$1"
    local git_root="$2"
    local build_dir="$3"

    # check that $commit is a commit hash, not something else
    local hash_="`git -C "$git_root" rev-parse "$commit"`"
    [ "x$commit" = "x$hash_" ]

    [ -e "$build_dir" ] || mkdir -p "`dirname "$build_dir"`"

    local logfile="${build_dir%/}-ogs-perf-build.log"

    if (
        set -x
        git -C "$git_root" switch -q --detach "$commit"
        # TODO always remove build dir to more easily recover from a previous
        # build failure?
        [ -e "$build_dir" ] || ogs_perf_init_build_dir "$git_root" "$build_dir"
        cd "$build_dir"
        ninja -j0 -k0
    ) &>"$logfile"
    then :
    else
        stat=$?
        error "tail of $logfile:"
        tail "$logfile"
        die "setup failed with status = $?."
    fi
}

# Sleeps conditionally for a dynamic amount of time.
#
# Arguments:
# - duration - the duration of a preceding task
# - max_no_sleep - the maximum duration (in s) such that we won't sleep
# - max_sleep - the maximum number of seconds to sleep
#
# The actual sleep time will be determined based on the passed arguments.
# The idea is that before starting the performance measurement a highly loaded
# system (e.g. from the build process) has some time to settle.
#
ogs_perf_conditional_sleep()
{
    local duration="$1"
    local max_no_sleep="$2"
    local max_sleep="$3"

    if [ "$duration" -gt "$max_no_sleep" ]
    then
        ((duration = duration * 2 / max_no_sleep))
        duration="$((duration < max_sleep ? duration : max_sleep))"
        echo "waiting for $duration s before continuing"
        sleep "$duration"
    fi
}

# Replaces certain characters in a git revision.
#
# Letters, numbers and -._ are kept. Any other character is replaced with an
# underscore.
#
# The resulting string can be used as a "nice" part of path names, as opposed to
# e.g. origin/master, HEAD^, my-remote/my-branch~5
#
ogs_perf_normalize_revision()
{
    local rev="$1"

    echo -n "$rev" | tr -s -c -- '-._[:alnum:]' '_'
}
