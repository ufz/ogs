#!/usr/bin/env bash

# Fail on added VTU files written in ASCII; existing ones are left alone.
if [[ -n ${PRE_COMMIT_FROM_REF-} ]]; then
    range=("$PRE_COMMIT_FROM_REF...$PRE_COMMIT_TO_REF")
else
    range=(--cached)
fi
files=$(git diff "${range[@]}" --diff-filter=A --name-only -z -- "$@" |
    xargs -0 -r grep -l 'format="ascii"' --)

if [[ -n $files ]]; then
    echo "Added VTU files in ASCII format; write them compressed binary:"
    echo "$files"
    exit 1
fi
