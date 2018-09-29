#!/bin/bash
if [[ $OSTYPE == darwin* ]]; then
    find ./Tests/Data -not -name "*.md5-stamp" -not -name "*expected*" -not -type d -exec bash -c 'mkdir -p ./tmp/`dirname {}`' \; -exec bash -c 'cp `grealpath {}` ./tmp/{}' \;
else
    find ./Tests/Data -not -name "*.md5-stamp" -not -name "*expected*" -not -type d -exec bash -c 'mkdir -p ./tmp/`dirname {}`' \; -exec bash -c 'cp `realpath {}` ./tmp/{}' \;
fi
pushd . > /dev/null
cd ./tmp/Tests/Data
tar zcf ogs6-data.tar.gz ./*
mv ogs6-data.tar.gz ../../../
zip -r -q ogs6-data .
mv ogs6-data.zip ../../../
popd > /dev/null
rm -rf tmp/
