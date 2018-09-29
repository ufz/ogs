#!/usr/bin/env bash

# A quick script to generate source tarball for the ogs edition.

usage="$0 <ParaView source dir>"
# A handy function for handling errors.
die () {
    echo >&2 "$@"
    exit 1
}

pv_dir="$1"

cp * $pv_dir/Catalyst
cd $pv_dir/Catalyst
./generate-tarballs.sh ogs
