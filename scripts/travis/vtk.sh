#!/bin/sh
set -e
if [ ! -d "$HOME/VTK-Install/include" ]; then
    cd $HOME
    wget http://ogsstorage.blob.core.windows.net/jenkins/ogs6-lib-sources/vtk-6.3.0-binary.tar.gz
    tar -xzf vtk-6.3.0.tar.gz
else
  echo 'Using cached vtk directory.';
fi
