#!/usr/bin/env bash

# set -e

# Check for uppercase letters in filename extensions
pattern=".*\.\w*([A-Z]+)\w*$"
files=$(echo $@ | grep -E "$pattern")

if [[ -n $files ]];
then
    echo "Attention!"
    echo "----------"
    echo "Found files that contain capital letters in the file extension."
    echo "Please rename the following files and commit again:"

    while read -r file; do
        echo -e '\E[0;32m'"$file"'\033[0m'
    done <<< "$files"
    # Abort commit
    exit 1
fi

exit 0
