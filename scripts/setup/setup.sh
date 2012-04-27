#!/usr/bin/env bash

# Parse options
while getopts ":a:d:s" opt; do
	case $opt in
		a)
			if [ "$OPTARG" == "x32" ]; then
				ARCHITECTURE="x32"
				WIN_ARCHITECTURE="x86"
			elif [ "$OPTARG" == "x64" ]; then
				ARCHITECTURE="x64"
				WIN_ARCHITECTURE="x64"
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
		s)
			echo "Qt SQL bindings enabled."
			QT_SQL=true
			;;
	esac
done

SOURCE_LOCATION=`pwd`
SOURCE_LOCATION="$SOURCE_LOCATION/../.."

source setup_prerequisites.sh

source setup_libraries.sh

source setup_builds.sh

if [ $QT_WAS_BUILT ]; then
	echo "Important note: Make sure to add the Qt environment variables!"
	if [ $QT_SQL ]; then
		echo "Important note: Make sure to add the instantclient directory to the PATH!"
	fi
fi