#!/bin/bash

# This script creates the quality assurance page and augments the
# input file parameter documentation by information about end-to-end
# tests in which the respective parameters are used.

if [ $# -ne 3 ]; then
    echo "USAGE: $0 SRCDIR BUILDDIR DATADIR" >&2
    exit 1
fi

if [ -d ".venv" ]; then
    source .venv/bin/activate
fi

set -e

srcdir="$1"
builddir="$2"
datadir="$3"

docauxdir="$builddir/DocAux"
doxdir="$docauxdir/dox"
toolsdir="$srcdir/scripts/doc"

param_cache="$docauxdir/documented-parameters-cache.txt"

qafile="$doxdir/project-file-doc-qa.dox"
check_quality_script="$toolsdir/check-project-params.py"

mkdir -p "$doxdir"

# Gather information about documented parameters.
"$toolsdir/get-project-params.sh" "$srcdir" \
    | "$toolsdir/normalize-param-cache.py" >"$param_cache"

# Document ctest project files
# and find out which tags and attributes are tested in which prj file and what
# is not tested at all.
"$toolsdir/linked-xml-file.py" "$datadir" "$docauxdir"

# Write QA information.
cat <<"EOF" >"$qafile"
/*! \page project_file_doc_qa OGS Input File Parameters&mdash;Quality Assurance

This page lists issues with the OGS input file parameter documentation.
If it is empty, there are no issues detected.

EOF

"$check_quality_script" "$docauxdir" "$srcdir" >>"$qafile" || true

cat <<EOF >>"$qafile"

*/
EOF

"$toolsdir/extract-media-properties-from-ctests.py" "$datadir" "$docauxdir"

# Finish parameter documentation dox files by appending auxiliary information,
# e.g., associated ctests, data type, etc.
"$toolsdir/append-xml-tags.py" prj "$datadir" "$docauxdir"
