#!/bin/sh
set -e
# check to see if cmake folder is empty
if [ ! -d "$HOME/eigen-eigen-dc6cfdf9bcec/Eigen" ]; then
    ZIP="3.2.9.zip"
    cd $HOME
    curl -L -o $ZIP http://bitbucket.org/eigen/eigen/get/$ZIP;
    unzip $ZIP;
else
    echo 'Using cached eigen directory.';
fi
