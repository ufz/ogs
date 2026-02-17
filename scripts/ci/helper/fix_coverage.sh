#!/usr/bin/env bash

SRC_ROOT="$1"

gcovr_args=""

if [[ -f "pytest-coverage.xml" ]]; then
  xmlstarlet ed \
    -u '//class/@filename' \
    -x 'concat("Applications/Python/", substring-after(., "site-packages/"))' \
    -u '//package/@name' -x 'substring-after(., "site-packages.")' \
    -d '//sources/source' \
    -s '//sources' -t elem -n source -v "$SRC_ROOT" \
    pytest-coverage.xml > pytest-coverage-rel.xml
  gcovr_args="--cobertura-add-tracefile pytest-coverage-rel.xml"
fi

if [[ -f "notebook-coverage.xml" ]]; then
  xmlstarlet ed \
    -u "//class/@filename[starts-with(., \"$SRC_ROOT/\")]" \
    -x "substring-after(., \"$SRC_ROOT/\")" \
    -u '//package/@name' -x 'substring-after(., "ogs.")' \
    -d '//sources/source' \
    -s '//sources' -t elem -n source -v "$SRC_ROOT" \
    notebook-coverage.xml > notebook-coverage-rel.xml
  gcovr_args="${gcovr_args} --cobertura-add-tracefile notebook-coverage-rel.xml"
fi

mkdir -p python-html
uvx gcovr \
  ${gcovr_args} \
  --html-details python-html/coverage.html \
  --exclude '.*/ogs/config.py' \
  -r ${SRC_ROOT} \
  -s