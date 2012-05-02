if [ -z $BUILD_LOCATION ]; then
	echo "You must specify a build directory (relative to the source directory)."
	echo "Aborting..."
	echo "Usage: $0 -d /build/directory"
	exit 1
fi

# Cleanup
rm -rf $BUILD_LOCATION
mkdir -p $BUILD_LOCATION && cd $BUILD_LOCATION