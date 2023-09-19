set -e

# Runs cppcheck with GitLab CI (CodeClimate) output
OUTPUT_FILE=${PROJECT_BINARY_DIR}/cppcheck.json
${CPPCHECK_TOOL_PATH} \
        --project=${PROJECT_BINARY_DIR}/compile_commands.json \
        --language=c++ \
        --std=c++${CMAKE_CXX_STANDARD} \
        --enable=all \
        --inconclusive \
        -j ${CPPCHECK_PARALLEL} \
        --inline-suppr \
        -i ${_cpm_dir} --suppress=*:${_cpm_dir}/* \
        -i ${PROJECT_SOURCE_DIR}/ThirdParty --suppress=*:${PROJECT_SOURCE_DIR}/ThirdParty* \
        -i ${PROJECT_SOURCE_DIR}/Tests --suppress=*:*Tests/* \
        --suppress=missingIncludeSystem \
        --template='{\n  "description": "{message}",\n  "severity": "info",\n    "location": {\n    "path": "{file}",\n    "lines": {\n      "begin": {line}\n    }\n  }\n},' \
        --output-file="$OUTPUT_FILE.tmp"

cat <<EOF >"$OUTPUT_FILE"
[
$(
  cat "$OUTPUT_FILE.tmp" | \
  `: strip source code absolute path` \
  sed 's|${PROJECT_SOURCE_DIR}/||' | \
  `: escape strings` \
  sed 's/string literal "\(.*\)" to/string literal \\"\1\\" to/g'
)
]
EOF

cat <<EOF >"$OUTPUT_FILE"
$(
  cat "$OUTPUT_FILE" | \
  `: remove last comma` \
  sed ${_last_sed_expression}
)
EOF

rm "$OUTPUT_FILE.tmp"

if [ -f ${Python_EXECUTABLE} ]; then
   ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/scripts/test/cppcheck_gen_hashes.py "$OUTPUT_FILE"
fi
