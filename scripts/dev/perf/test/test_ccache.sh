#!/bin/sh

# These tests check ccache effectiveness under certain conditions.

set -eu
set -o pipefail
shopt -s inherit_errexit

# make set -x more useful, cf. https://web.archive.org/web/20230401201759/https://wiki.bash-hackers.org/scripting/debuggingtips
export PS4='+<(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}()> }'

exists() { [ -e "$1" ]; }

### some definitions
#
today="`date +%F`"

# to make dirs/branches unique
uuid="`uuidgen`"
#
# temporary branch
branch="test-master-$today-$uuid"

worktree="debug-tmpworktree-`date +%F_%H%M%S`_$uuid"

counter=0

generate_source_and_build()
{
    local worktree="`realpath "$1"`"
    shift
    local srcdir="$worktree/src"
    local builddir="$worktree/build"

    mkdir -p "$worktree"
    mkdir "$srcdir"
    mkdir "$builddir"

    cat >"$srcdir/test.cpp" <<EOF
#include <iostream>

int main()
{
    std::cout << "$uuid" << '\n';
}
EOF

    export CCACHE_BASEDIR="$worktree"
    export CCACHE_SLOPPINESS=pch_defines,time_macros,mtime
    export CCACHE_NOHASHDIR=''
    export CCACHE_STATSLOG="${builddir%/}_ccache_stats_`date +%F`.log"
    export CCACHE_DEBUG=1

    export CXXFLAGS="-Wp,-D_GLIBCXX_ASSERTIONS -fno-omit-frame-pointer -mno-omit-leaf-frame-pointer"

    (
        cd "$builddir"  # is important for cache hits!
        ccache "$CXX" -O3 -Wall \
            "$@" \
            -c "$srcdir/test.cpp" -o "$builddir/test.cpp.o"

        statslog="$(ccache --show-log-stats)"
        echo "$statslog"
    )
}

get_percent_hits()
{
    echo "$1" | sed -n -e "/Hits:/{s/^ *Hits:.*( *\([0-9]*\)\..*%)/\1/p; q}"
}

success=1


echo "==============================="
ccache --version
echo "==============================="
g++ --version
echo "==============================="
clang++ --version
echo "==============================="

# set -x

for CXX in g++ clang++
do
    export CXX="$CXX"

    for MARCH in x86-64 native
    do
        # Skip clang++ with native, Remove when ccache is released with
        # https://github.com/ccache/ccache/issues/1658.
        [ "$CXX" = "clang++" ] && [ "$MARCH" = "native" ] && continue

        MARCH="-march=$MARCH"

        for debug in "" -g
        do
            ((++counter))

            # testee - first build, cache miss expected
            statslog="`generate_source_and_build "$worktree-$counter" "$MARCH" $debug`"

            perc_hits="`get_percent_hits "$statslog"`"
            echo "$CXX [$MARCH] [$debug] $perc_hits% cache hits"

            if [ 0 -ne "$perc_hits" ]; then
                echo "  FAIL"
                success=0
            fi

            ((++counter))

            # testee - second build in new dir, cache hit expected
            statslog="`generate_source_and_build "$worktree-$counter" "$MARCH" $debug`"

            perc_hits="`get_percent_hits "$statslog"`"
            echo "$CXX [$MARCH] [$debug] $perc_hits% cache hits"

            if [ 100 -ne "$perc_hits" ]; then
                echo "  FAIL"
                success=0
            fi
        done
    done
done

[ 1 -eq "$success" ]
