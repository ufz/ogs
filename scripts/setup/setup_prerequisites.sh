#!/usr/bin/env bash

# Check for prerequisites
CMAKE_LOCATION=`which cmake`
if [ -z "$CMAKE_LOCATION" ]; then
	echo "CMake not found! Aborting..."
	exit 1
fi

source $SOURCE_LOCATION/scripts/base/download_file_with_cmake.sh
source $SOURCE_LOCATION/scripts/base/configure_compiler.sh

## Windows specific
if [ "$OSTYPE" == 'msys' ]; then

	mkdir -vp ~/bin

	# 7-zip
	SEVENZIP_LOCATION=`which 7za`
	if [ ! -z "$SEVENZIP_LOCATION" ]; then
		echo "7-zip found."
	else
		cd ~/bin
		download_file http://dl.dropbox.com/u/5581063/7za.exe ./7za.exe e92604e043f51c604b6d1ac3bcd3a202
	fi
	SEVENZIP_LOCATION=`which 7za`
	if [ -z "$SEVENZIP_LOCATION" ]; then
		echo "7-zip not downloaded! Aborting..."
		exit 1
	fi

	# jom
	JOM_LOCATION=`which jom`
	if [ ! -z "$JOM_LOCATION" ]; then
		echo "jom found."
	else
		cd ~/bin
		download_file http://dl.dropbox.com/u/5581063/jom.exe ./jom.exe 335428f223d36f0a39faaa845222346d
	fi
	JOM_LOCATION=`which jom`
	if [ -z "$JOM_LOCATION" ]; then
		echo "jom not downloaded! Aborting..."
		exit 1
	fi

fi


cd "$SOURCE_LOCATION/scripts/setup"