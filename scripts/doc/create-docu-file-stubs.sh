#!/bin/bash

# This script creates missing file in the OGS input file documentation source tree.
# It expects input from get-project-params.sh
#
# The working directory of this script must be the root of the OGS sources.
#
# Example:
# scripts/doc/get-project-params.sh . | scripts/doc/create-docu-file-stubs.sh

base="Documentation/ProjectFile"

while IFS=":" read -r fn lno content; do
    [ "$content" = "${content#*//!}" ] && continue
    tag_name="$(echo "$content" \
        | sed -n -e 'sX^\s*//! \\ogs_file_\(param\|attr\)\(_special\)\?{\([A-Za-z_0-9]\+\)}$X\1 \3Xp')"
    [ -z "$tag_name" ] && continue
    param_or_attr="${tag_name%% *}"
    tag_name="${tag_name#* }"
    tag_name="${tag_name//__/\/}"
    echo "$param_or_attr $base/$tag_name"
done \
| sort -r \
| while read param_or_attr path; do
    dn="`dirname "$path"`"
    bn="`basename "$path"`"
    # echo "$param_or_attr $path"

    if [ ! -d "$dn" ]; then
        mkdir -p "$dn"

        bdn="`basename "$dn"`"
        if [ "`expr match "$bdn" '^[A-Z]'`" -eq 0 ] && [ ! -f "$dn/i_$bdn.md" ]; then
            echo "creating $dn/i_$bdn.md"
            echo '\ogs_missing_documentation' >"$dn/i_$bdn.md"
        elif [ "`expr match "$bdn" '^[A-Z]'`" -ne 0 ] && [ ! -f "$dn/c_$bdn.md" ]; then
            echo "creating $dn/c_$bdn.md"
            echo '\ogs_missing_documentation' >"$dn/c_$bdn.md"
        fi
    fi

    if [ -d "$path" ]; then
        if [ "`expr match "$bn" '^[A-Z]'`" -eq 0 ] && [ ! -f "$path/i_$bn.md" ]; then
            echo "creating $path/i_$bn.md"
            echo '\ogs_missing_documentation' >"$path/i_$bn.md"
        elif [ "`expr match "$bn" '^[A-Z]'`" -ne 0 ] && [ ! -f "$path/c_$bn.md" ]; then
            echo "creating $path/c_$bn.md"
            echo '\ogs_missing_documentation' >"$path/c_$bn.md"
        fi
    elif [ "$param_or_attr" = param ] && [ ! -f "$dn/t_$bn.md" ]; then
        echo "creating $dn/t_$bn.md"
        echo '\ogs_missing_documentation' >"$dn/t_$bn.md"
    elif [ "$param_or_attr" = attr ] && [ ! -f "$dn/a_$bn.md" ]; then
        echo "creating $dn/a_$bn.md"
        echo '\ogs_missing_documentation' >"$dn/a_$bn.md"
    # else
    #     echo "OK $path"
    fi

    # if [ -d "$path" ] && [ -f "$path.md" ]; then
    #     echo "ERROR: both $path and $path.md exist!" >&2
    # fi
done
