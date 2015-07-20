#!/bin/sh
set -e
# check to see if cmake folder is empty
if [ ! -d "$HOME/cmake-3.1.1-Linux-x86_64/bin" ]; then
	cd $HOME
	wget http://www.cmake.org/files/v3.1/cmake-3.1.1-Linux-x86_64.tar.gz;
	tar -xzvf cmake-3.1.1-Linux-x86_64.tar.gz;
else
	echo 'Using cached cmake directory.';
fi
