#!/usr/bin/env bash

# Parse options
while getopts ":a:d" opt; do
	case $opt in
		a)
			if [ "$OPTARG" == "x32" ]; then
				ARCHITECTURE="x32"
			elif [ "$OPTARG" == "x64" ]; then
				ARCHITECTURE="x64"
			else
				echo "$OPTARG is not a valid argument. Specify x32 or x64."
				exit 1
			fi
			;;
		\?)
			echo "Invalid option: -$OPTARG"
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument."
			exit 1
			;;
		d)
			echo "Third party library debug builds enabled."
			LIB_DEBUG=true
			;;
	esac
done

SOURCE_LOCATION=`pwd`
SOURCE_LOCATION="$SOURCE_LOCATION/../.."

source setup_prerequisites.sh

source setup_libraries.sh