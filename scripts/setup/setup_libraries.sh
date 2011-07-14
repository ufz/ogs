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
			" > build_qt.bat

			$COMSPEC \/k build_qt.bat
		fi
		
		export PATH=$PATH:$SOURCE_LOCATION/../libs/$QT_VERSION/bin
		
	else
		echo "Qt already installed in $QMAKE_LOCATION"
	fi
	
	cd "$SOURCE_LOCATION/../libs"
	
	# Install VTK
	#http://www.vtk.org/files/release/5.6/vtk-5.6.1.tar.gz
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
fi