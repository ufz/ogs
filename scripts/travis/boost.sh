#!/bin/sh
set -e
# check to see if boost folder is empty
if [ ! -d "$HOME/boost_1_56_0/boost" ]; then
	TAR="boost_1_56_0.tar.gz"
	cd $HOME
	curl -L -o $TAR https://sourceforge.net/projects/boost/files/boost/1.56.0/$TAR/download;
	tar -xzf $TAR;
else
	echo 'Using cached boost directory.';
fi
