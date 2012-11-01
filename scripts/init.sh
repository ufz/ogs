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
		download_file http://sourceforge.net/projects/uncrustify/files/uncrustify/uncrustify-0.59/uncrustify-0.59-win32.zip/download uncrustify-0.59-win32.zip
		7za x uncrustify-0.59-win32.zip
		mv uncrustify-0.59-win32/uncrustify.exe ~/bin
		rm -r uncrustify-0.59-win32*
		if [ ! -f ~/bin/uncrustify.exe ]; then
			echo "Error downloading uncrustify! Git hooks not set."
		fi
	else
		echo "Please install uncrustify (http://uncrustify.sourceforge.net) to setup git hooks."
		exit 1
	fi
fi

cd $DIR
# Enable git hooks
cd .git/hooks
git init
cd ../..
git fetch origin
cd .git/hooks
git pull .. remotes/origin/hooks
cd $DIR

# Set git configs for running uncrustify
git config --bool hooks.uncrustify
git config hooks.uncrustify.conf $SCRIPT_DIR/style/uncrustify.cfg
