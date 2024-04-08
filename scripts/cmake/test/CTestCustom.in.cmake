file(STRINGS "@PROJECT_BINARY_DIR@/CTestTestfile.cmake" LINES)

# overwrite the file....
file(WRITE "@PROJECT_BINARY_DIR@/CTestTestfile.cmake" "")

# loop through the lines,
foreach(line IN LISTS LINES)
    # remove unwanted parts
    string(REGEX REPLACE ".*_deps/.*" "" STRIPPED "${line}")
    # and write the (changed) line ...
    file(APPEND "@PROJECT_BINARY_DIR@/CTestTestfile.cmake" "${STRIPPED}\n")
endforeach()
