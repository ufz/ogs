# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DIR=$SCRIPT_DIR/..
source $SCRIPT_DIR/base/download_file_with_cmake.sh

### Init git hooks ###

# Check for uncrustify
UNCRUSTIFY_LOCATION=`which uncrustify`
if [ -z "$UNCRUSTIFY_LOCATION" ]; then
	if [ "$OSTYPE" == 'msys' ]; then
		mkdir -vp ~/bin
		export PATH=$PATH:~/bin
		download_file http://sourceforge.net/projects/uncrustify/files/uncrustify/uncrustify-0.59/uncrustify-0.59-win32.zip/download ./uncrustify-0.59-win32.zip
		7za x uncrustify-0.59-win32.zip
		mv uncrustify-0.59-win32/uncrustify.exe ~/bin
		rm -r uncrustify-0.59-win32*
		if [ ! -f ~/bin/uncrustify.exe ]; then
			echo "Error downloading uncrustify! Git hooks not set."
		else
			UNCRUSTIFY_LOCATION=~/bin/uncrustify.exe
		fi
	else
		echo "Please install uncrustify (http://uncrustify.sourceforge.net) to setup git hooks."
		exit 1
	fi
fi

cd $DIR
# Enable git hooks
# Does not work on Windows at the moment
if [ ! "$OSTYPE" == 'msys' ]; then
	cd .git/hooks
	git init
	cd ../..
	git fetch origin
	cd .git/hooks
	git pull .. remotes/origin/hooks
	cd $DIR

	# Set git configs for running uncrustify
	git config --add --bool hooks.uncrustify true
	git config hooks.uncrustify.path $UNCRUSTIFY_LOCATION
	git config hooks.uncrustify.conf $SCRIPT_DIR/style/uncrustify.cfg
fi
