# Check Visual Studio version and setup CMake generator
if [ -z "$VS110COMNTOOLS" ]; then
	if [ -z "$VS100COMNTOOLS" ]; then
		if [ -z "$VS90COMNTOOLS" ]; then
			if [ -z "$VS80COMNTOOLS" ]; then
				echo "Error: Visual Studio not found"
				exit 1
			else
				WIN_DEVENV_PATH="$VS80COMNTOOLS..\\IDE"
				CMAKE_GENERATOR="Visual Studio 8 2005"
			fi
		else
			WIN_DEVENV_PATH="$VS90COMNTOOLS..\\IDE"
			CMAKE_GENERATOR="Visual Studio 9 2008"
		fi
	else
		WIN_DEVENV_PATH="$VS100COMNTOOLS..\\IDE\\"
		CMAKE_GENERATOR="Visual Studio 10"
	fi
else
	WIN_DEVENV_PATH="$VS110COMNTOOLS..\\IDE\\"
	CMAKE_GENERATOR="Visual Studio 11"
	if [ "$ARCHITECTURE" == "x64" ]; then
		WIN_ARCHITECTURE="x86_amd64"
	fi
fi

if [ "$ARCHITECTURE" == "x64" ]; then
	CMAKE_GENERATOR="$CMAKE_GENERATOR Win64"
fi

DEVENV_EXE="${WIN_DEVENV_PATH}devenv"

# Replace backslashes in WIN_DEVENV_PATH
DEVENV_PATH=$(echo "$WIN_DEVENV_PATH" | awk '{ gsub(/\\/, "/"); print }')
DEVENV_PATH=$(echo "$DEVENV_PATH" | awk '{ gsub(/C:\//, "/c/"); print }')

echo "Visual Studio found: $DEVENV_PATH"
echo "Devenv: $DEVENV_EXE"
echo "CMake Generator: $CMAKE_GENERATOR"
export PATH=$PATH:$DEVENV_PATH