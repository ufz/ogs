#!/bin/sh
set -e
if [ ! -d "$HOME/VTK-Install/include" ]; then
	cd $HOME
	wget --no-check-certificate http://opengeosys.s3.amazonaws.com/ogs6-lib-sources/vtk-6.3.0.tar.gz
	tar -xzf vtk-6.3.0.tar.gz
else
  echo 'Using cached vtk directory.';
fi
