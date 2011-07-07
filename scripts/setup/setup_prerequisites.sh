#!/usr/bin/env bash

## Windows specific
if [ "$OSTYPE" == 'msys' ]; then
	
	# Check Visual Studio version
	if [ -z "$VS80COMNTOOLS" ]; then
		if [ -z "$VS90COMNTOOLS" ]; then
			if [ -z "VS100COMNTOOLS" ]; then
				echo "Error: Visual Studio not found"
			else
				DEVENV="$VS100COMNTOOLS..\IDE\devenv"
			fi
		else
			DEVENV="$VS90COMNTOOLS..\IDE\devenv"
		fi
	else
		DEVENV="$VS80COMNTOOLS..\IDE\devenv"
	fi

	echo "Visual Studio found: $DEVENV"
fi
