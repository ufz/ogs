#!/usr/bin/env bash

EXE=$(command -v cmake-format)
if [ -z "$EXE" ]; then
    echo "cmake-format not found; skipping check."
    exit 0
fi


MOD_FILES=""

MOD_FILES=()
for FILE in "$@"
do
    cmake-format -i "$FILE"
    MODIFIED=$(git status --porcelain "$FILE" | head -c 2 | tail -c 1)
    # echo $FILE: $MODIFIED
    if [ "$MODIFIED" = "M" ]; then
        MOD_FILES+=("$FILE")
        echo "Fixed $FILE"
    fi
done

if [ -z "$MOD_FILES" ]; then
    exit 0
fi

echo "CMake files have been modified."
echo "Add them to the commit!"
exit 1
