#!/usr/bin/env bash

## Windows specific
if [ "$OSTYPE" == 'msys' ]; then
	
	# Check Visual Studio version and setup CMake generator
	if [ -z "$VS80COMNTOOLS" ]; then
		if [ -z "$VS90COMNTOOLS" ]; then
			if [ -z "VS100COMNTOOLS" ]; then
				echo "Error: Visual Studio not found"
				exit 1
			else
				WIN_DEVENV_PATH="$VS100COMNTOOLS..\\IDE\\"
				CMAKE_GENERATOR="Visual Studio 10"
			fi
		else
			WIN_DEVENV_PATH="$VS90COMNTOOLS..\\IDE"
			CMAKE_GENERATOR="Visual Studio 9 2008"
		fi
	else
		WIN_DEVENV_PATH="$VS80COMNTOOLS..\\IDE"
		CMAKE_GENERATOR="Visual Studio 8 2005"
	fi
	
	if [ "$ARCHITECTURE" == "x64" ]; then
		CMAKE_GENERATOR="$CMAKE_GENERATOR Win64"
	fi
	
	# Replace backslashes in WIN_DEVENV_PATH
	DEVENV_PATH=$(echo "$WIN_DEVENV_PATH" | awk '{ gsub(/\\/, "/"); print }')
	DEVENV_PATH=$(echo "$DEVENV_PATH" | awk '{ gsub(/C:\//, "/c/"); print }')
	
	echo "Visual Studio found: $DEVENV_PATH"
	echo "CMake Generator: $CMAKE_GENERATOR"
	export PATH=$PATH:$DEVENV_PATH
	
	# 7-zip
	SEVENZIP_LOCATION=`which 7za`
	if [ ! -z "$SEVENZIP_LOCATION" ]; then
		echo "7-zip found."
	else
		cd ~/bin
		wget --no-check-certificate https://github.com/downloads/ufz/devguide/7za.exe
	fi

fi

cd "$SOURCE_LOCATION/scripts/setup"