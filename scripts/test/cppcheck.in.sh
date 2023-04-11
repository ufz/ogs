set -e

# Runs cppcheck with GitLab CI (CodeClimate) output
OUTPUT_FILE=${PROJECT_BINARY_DIR}/cppcheck.json
${CPPCHECK_TOOL_PATH} \
        --project=${PROJECT_BINARY_DIR}/compile_commands.json \
        --language=c++ \
        --std=c++20 \
        --enable=all \
        --inconclusive \
        -j ${CPPCHECK_PARALLEL} \
        ${_cpp_check_ingore} \
        --inline-suppr \
        --suppress=*:*/usr/local\* \
        --suppress=*:*cpm\* \
        --suppress=*:*Tests\* \
        --template='{\n  "description": "{message}",\n  "severity": "info",\n    "location": {\n    "path": "{file}",\n    "lines": {\n      "begin": {line}\n    }\n  }\n},' \
        --output-file="$OUTPUT_FILE.tmp"

cat <<EOF >"$OUTPUT_FILE"
[
$(
  cat "$OUTPUT_FILE.tmp" | \
  `: strip source code absolute path` \
  sed 's|${PROJECT_SOURCE_DIR}/||' | \
  `: escape strings` \
  sed 's/string literal "\(.*\)" to/string literal \\"\1\\" to/g' | \
  `: remove last comma` \
  sed '$s/,$//'
)
]
EOF

rm "$OUTPUT_FILE.tmp"

if [ -f ${Python_EXECUTABLE} ]; then
    ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/scripts/test/cppcheck_gen_hashes.py "$OUTPUT_FILE"
fi
