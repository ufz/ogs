#!/usr/bin/env bash

cd "$SOURCE_LOCATION/../"
mkdir -vp libs
cd libs

QMAKE_LOCATION=`which qmake`

## Windows specific
if [ "$OSTYPE" == 'msys' ]; then
	if [ -z "$QMAKE_LOCATION" ]; then
		# Install Qt
		QT_VERSION="qt-everywhere-opensource-src-4.7.3"
		
		if [ ! -d $QT_VERSION ]; then
			# Download and extract
			wget http://get.qt.nokia.com/qt/source/$QT_VERSION.zip -O ./$QT_VERSION.zip
			7za x $QT_VERSION.zip
			rm $QT_VERSION.zip
		elif [ -f $QT_VERSION/bin/qmake.exe -a -f $QT_VERSION/bin/QtGui4.dll ]; then
			# Already installed
			QT_FOUND=true
		fi
		
		if [ $QT_FOUND ]; then
			echo "Qt already installed in $SOURCE_LOCATION/../$QT_VERSION"
		else
			# Compile
			
			# TODO: -mp flag for multiprocessor compiling?
			if [ $LIB_DEBUG ]; then
				QT_CONFIGURATION="-debug-and-release"
			else
				QT_CONFIGURATION="-release"
			fi
			
			cd $QT_VERSION

			echo " \
			\"$WIN_DEVENV_PATH\\..\\..\\VC\\vcvarsall.bat\" $WIN_ARCHITECTURE &&\
			echo y | configure -opensource -nomake demos -nomake examples $QT_CONFIGURATION &&\
			nmake && nmake clean &&\
			exit\
			" > build.bat

			$COMSPEC \/k build.bat
		fi
		
		export PATH=$PATH:$SOURCE_LOCATION/../libs/$QT_VERSION/bin
		
	else
		echo "Qt already installed in $QMAKE_LOCATION"
	fi
	
	
	# Install VTK
	cd "$SOURCE_LOCATION/../libs"
	VTK_VERSION="vtk-5.6.1"
	if [ ! -d $VTK_VERSION ]; then
		# Download, extract, rename
		wget http://www.vtk.org/files/release/5.6/$VTK_VERSION.tar.gz -O ./$VTK_VERSION.tar.gz
		tar -xf $VTK_VERSION.tar.gz
		mv VTK/ $VTK_VERSION/
		rm $VTK_VERSION.tar.gz
	# Check for existing installation
	elif [ -f $VTK_VERSION/build/bin/Release/QVTK.lib -a -f $VTK_VERSION/build/bin/Release/vtkRendering.lib ]; then
		if [ $LIB_DEBUG ]; then
			if [ -f $VTK_VERSION/build/bin/Debug/QVTK.lib -a -f $VTK_VERSION/build/bin/Debug/vtkRendering.lib ]; then
				VTK_FOUND=true
			fi
		else
			VTK_FOUND=true
		fi
	fi
	
	if [ $VTK_FOUND ]; then
		echo "VTK already installed in $SOURCE_LOCATION/../$VTK_VERSION"
	else
		# Compile
		cd $VTK_VERSION
		mkdir -vp build
		cd build
		cmake .. -DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=OFF -DVTK_USE_GUISUPPORT=ON -DVTK_USE_QT=ON -DVTK_USE_QVTK_QTOPENGL=ON -G "$CMAKE_GENERATOR"
		cmake ..
		echo "PATH: $PATH"
		qmake -v
		$COMSPEC \/c "devenv.com VTK.sln /Build Release"
		$COMSPEC \/c "devenv VTK.sln /Build Release /Project QVTK"
		
		if [ $LIB_DEBUG ]; then
			$COMSPEC \/c "devenv VTK.sln /Build Debug"
			$COMSPEC \/c "devenv VTK.sln /Build Debug /Project QVTK"
		fi
	fi
	
	# Install shapelib
	cd "$SOURCE_LOCATION/../libs"
	SHAPELIB_VERSION="shapelib-1.3.0b2"
	if [ ! -d $SHAPELIB_VERSION ]; then
		# Download, extract
		wget http://download.osgeo.org/shapelib/$SHAPELIB_VERSION.tar.gz
		tar -xf $SHAPELIB_VERSION.tar.gz
		rm -rf $SHAPELIB_VERSION.tar.gz
	elif [ -f $SHAPELIB_VERSION/shapelib.lib ]; then
		SHAPELIB_FOUND=true
	fi
	
	if [ $SHAPELIB_FOUND ]; then
		echo "Shapelib already installed in $SOURCE_LOCATION/../$SHAPELIB_VERSION"
	else
		# Compile
		cd $SHAPELIB_VERSION
		
		echo " \
		\"$WIN_DEVENV_PATH\\..\\..\\VC\\vcvarsall.bat\" $WIN_ARCHITECTURE &&\
		nmake /f makefile.vc &&\
		exit\
		" > build.bat

		$COMSPEC \/k build.bat
	fi
	
	# Install libtiff
	cd "$SOURCE_LOCATION/../libs"
	LIBTIFF_VERSION="tiff-3.9.5"
	if [ ! -d $LIBTIFF_VERSION ]; then
		# Download, extract
		wget ftp://ftp.remotesensing.org/pub/libtiff/$LIBTIFF_VERSION.tar.gz
		tar -xf $LIBTIFF_VERSION.tar.gz
		rm -rf $LIBTIFF_VERSION.tar.gz
	elif [ -f $LIBTIFF_VERSION/libtiff/libtiff.lib ]; then
		LIBTIFF_FOUND=true
	fi
	
	if [ $LIBTIFF_FOUND ]; then
		echo "Libtiff already installed in $SOURCE_LOCATION/../$LIBTIFF_VERSION"
	else
		# Compile
		cd $LIBTIFF_VERSION
		
		echo " \
		\"$WIN_DEVENV_PATH\\..\\..\\VC\\vcvarsall.bat\" $WIN_ARCHITECTURE &&\
		nmake /f Makefile.vc lib &&\
		exit\
		" > build.bat

		$COMSPEC \/k build.bat
	fi
	
	# Install libgeotiff
	cd "$SOURCE_LOCATION/../libs"
	LIBGEOTIFF_VERSION="libgeotiff-1.3.0"
	if [ ! -d $LIBGEOTIFF_VERSION ]; then
		# Download, extract
		wget http://download.osgeo.org/geotiff/libgeotiff/$LIBGEOTIFF_VERSION.tar.gz
		tar -xf $LIBGEOTIFF_VERSION.tar.gz
		rm -rf $LIBGEOTIFF_VERSION.tar.gz
	elif [ -f $LIBGEOTIFF_VERSION/geotiff.lib ]; then
		LIBGEOTIFF_FOUND=true
	fi
	
	if [ $LIBGEOTIFF_FOUND ]; then
		echo "Libgeotiff already installed in $SOURCE_LOCATION/../$LIBGEOTIFF_VERSION"
	else
		# Compile
		cd $LIBGEOTIFF_VERSION
		
		# Download modified makefile
		if [ ! -f makefile_mod.vc ]; then
			wget --no-check-certificate https://gist.github.com/raw/1088657/0b846a9cdc529681bfb34be37dfba5d1a31dc419/makefile_mod.vc
		fi
		
		echo " \
		\"$WIN_DEVENV_PATH\\..\\..\\VC\\vcvarsall.bat\" $WIN_ARCHITECTURE &&\
		nmake /f makefile_mod.vc geotiff.lib&&\
		exit\
		" > build.bat

		$COMSPEC \/k build.bat
	fi
fi