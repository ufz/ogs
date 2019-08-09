#!/usr/bin/env bash

set -e

BINARY_FILES=""
LFS_FILES=""
CHANGED_FILES=$(git diff --cached --name-only --diff-filter=ACMRTXBU)

while read -r FILE; do
    LFS_FILE=$(git check-attr filter "$FILE" | grep 'filter: lfs$' | sed -e 's/: filter: lfs//')
    if [ ! -z "$LFS_FILE" ]; then
        LFS_FILES="$LFS_FILES $LFS_FILE"
    fi
done <<< "$CHANGED_FILES"

if [ -z "$LFS_FILES" ]; then
    exit 0
fi

while read -r FILE; do
    SOFT_SHA=$(git hash-object -w "$FILE")
    RAW_SHA=$(git hash-object -w --no-filters "$FILE")

    if [ $SOFT_SHA == $RAW_SHA ]; then
        BINARY_FILES="$FILE\n$BINARY_FILES"
    fi
done <<< "$LFS_FILES"

if [[ -n "$BINARY_FILES" ]]; then
    echo "Attention!"
    echo "----------"
    echo "You tried to commit binary files:"
    echo -e "\x1B[31m$BINARY_FILES\x1B[0m"
    echo "Make sure you have git-lfs installed properly!"
    echo "Revert your changes and commit those files with git-lfs!"
    echo "----------"
    exit 1
fi
