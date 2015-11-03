#!/bin/sh
set -e
if [ ! -d "$HOME/VTK-Install/include" ]; then
	cd $HOME
	wget --no-check-certificate http://www.opengeosys.org/images/dev/vtk-6.1.0.tar.gz
	tar -xzf vtk-6.1.0.tar.gz
else
  echo 'Using cached vtk directory.';
fi
