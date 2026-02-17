#!/usr/bin/env bash

SRC_ROOT="$1"

xmlstarlet ed \
  -u '//class/@filename' \
  -x 'concat("Applications/Python/", substring-after(., "site-packages/"))' \
  -u '//package/@name' -x 'substring-after(., "site-packages.")' \
  -d '//sources/source' \
  -s '//sources' -t elem -n source -v "$SRC_ROOT" \
  pytest-coverage.xml > pytest-coverage-rel.xml

xmlstarlet ed \
  -u '//class/@filename' \
  -x 'concat("Tests/Data/", substring-after(., "Tests/Data/"))' \
  -u '//package/@name' -x 'substring-after(., "ogs.")' \
  -d '//sources/source' \
  -s '//sources' -t elem -n source -v "$SRC_ROOT" \
  notebook-coverage.xml > notebook-coverage-rel.xml
