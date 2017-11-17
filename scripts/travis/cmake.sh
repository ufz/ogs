#!/bin/sh
set -e
# check to see if cmake folder is empty
if [ ! -d "$HOME/cmake-3.1.1-Linux-x86_64/bin" ]; then
    CMAKE_TAR="cmake-3.1.1-Linux-x86_64.tar.gz"
    cd $HOME
    curl -L -o $CMAKE_TAR http://www.cmake.org/files/v3.1/$CMAKE_TAR;
    tar -xzf $CMAKE_TAR;
else
    echo 'Using cached cmake directory.';
fi
