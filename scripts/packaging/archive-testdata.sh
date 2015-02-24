#!/bin/bash
if [[ $OSTYPE == darwin* ]]; then
	find ./Tests/Data -not -name "*.md5-stamp" -not -name "*expected*" -not -type d -exec mkdir -p ./tmp/{} \; -exec cp `grealpath {}` ./tmp/{} \;
else
	find ./Tests/Data -not -name "*.md5-stamp" -not -name "*expected*" -not -type d -exec mkdir -p ./tmp/{} \; -exec cp `realpath {}` ./tmp/{} \;
fi
pushd . > /dev/null
cd ./tmp/Tests/Data
tar zcf ogs6-data.tar.gz ./*
mv ogs6-data.tar.gz ../../../
popd > /dev/null
rm -rf tmp/
