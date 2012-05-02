NUM_PROCESSORS=""
if [ "$OSTYPE" == 'msys' ]; then
	source $SOURCE_LOCATION/scripts/base/configure_win_vs.sh
else
	CMAKE_GENERATOR="Unix Makefiles"
	if [[ "$OSTYPE" == darwin* ]]; then
		NUM_PROCESSORS=`sysctl hw.ncpu | awk '{print $2}'`
	else
		NUM_PROCESSORS=`cat /proc/cpuinfo | grep processor | wc -l`
	fi
fi