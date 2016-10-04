#!/bin/bash

# This script traverses the OGS input file documentation source tree
# making a list of all input file parameter related Doxygen commands,
# i.e., those beginning with \ogs, and of all ConfigTree related C++
# code.

if [ $# -ne 1 ]; then
    echo "USAGE: ${0##*/} SRCDIR" >&2
    exit 1
fi

srcdir="$1"

#color="--color=always"
color=""

cat <<"EOF" \
| grep -r "$srcdir" \
    --include '*.h' \
    --include '*.cpp' \
    --exclude-dir '.git' \
    --exclude-dir 'Tests' \
    --exclude 'ConfigTree*.*' \
    -f - -r -n -o $color \
| cut -c $((${#srcdir} + 2))-
//! \\ogs_file_\(param\|attr\){[A-Za-z_0-9]\+}\( \\todo .*\)\?$
//! \\ogs_file_special$
//! \\ogs_file_\(param\|attr\)_special{[A-Za-z_0-9]\+}\( \\todo .*\)\?$
checkConfigParameter[^)]*)\?
getConfigAttribute[^)]*)\?
getConfigParameter[^)]*)\?
getConfigSubtree[^)]*)\?
ignoreConfigAttribute[^)]*)\?
ignoreConfigParameter[^)]*)\?
peekConfigParameter[^)]*)\?
EOF

# format as table:
# | sed -e 's_::_@@_g' -e's_:\s\+_:_' | column -t -s: | sed -e 's_@@_::_g'
