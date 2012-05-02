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