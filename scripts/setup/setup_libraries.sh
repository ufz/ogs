#!/usr/bin/env bash

LIBS_LOCATION="$SOURCE_LOCATION/../libs"

mkdir -vp $LIBS_LOCATION
cd $LIBS_LOCATION
mkdir -vp include
mkdir -vp lib

QMAKE_LOCATION=`which qmake`

QT_VERSION="qt-everywhere-opensource-src-4.8.3"
VTK_VERSION="vtk-5.10.0"
SHAPELIB_VERSION="shapelib-1.3.0"
LIBGEOTIFF_VERSION="libgeotiff-1.3.0"
INSTANTCLIENT_VERSION="instantclient_11_2"
METIS_VERSION="metis-5.0.2"
BOOST_VERSION="1.53.0"
BOOST_VERSION_UNDERSCORE=${BOOST_VERSION//./_}

## Windows specific
if [ "$OSTYPE" == 'msys' ]; then
	if [ -z "$QMAKE_LOCATION" ]; then
		# Install Qt
		if [ ! -d qt ]; then
			# Download and extract
			download_file http://releases.qt-project.org/qt4/source/$QT_VERSION.zip ./$QT_VERSION.zip
			7za x $QT_VERSION.zip
			mv $QT_VERSION/ qt/
			rm $QT_VERSION.zip

		elif [ -f qt/bin/qmake.exe -a -f qt/bin/QtGui4.dll -a -f qt/bin/QtGuid4.dll ]; then
			# Already installed
			QT_FOUND=true
		fi

		if [ $QT_FOUND ]; then
			echo "Qt already installed in $LIBS_LOCATION/qt"
		else
			# Compile
			QT_CONFIGURATION="-debug-and-release"

			# Get instantclient
			QT_SQL_ARGS=""
			if [ $QT_SQL ]; then
				if [ ! -d instantclient ]; then
					if [ "$ARCHITECTURE" == "x64" ]; then
						download_file http://dl.dropbox.com/u/5581063/instantclient_11_2_x64.zip ./instantclient_11_2_x64.zip 015bd1b163571988cacf70e7d6185cb5
						7za x instantclient_11_2_x64.zip
						mv instantclient_11_2/ instantclient/
						rm instantclient_11_2_x64.zip
					fi
				fi
				QT_SQL_ARGS="-qt-sql-oci -I %cd%\..\instantclient\sdk\include -L %cd%\..\instantclient\sdk\lib\msvc"
			fi

			cd qt

			echo " \
			\"$WIN_DEVENV_PATH\\..\\..\\VC\\vcvarsall.bat\" $WIN_ARCHITECTURE &&\
			echo y | configure -opensource -no-accessibility -no-dsp -no-vcproj -no-phonon -no-webkit -no-scripttools -nomake demos -nomake examples $QT_CONFIGURATION $QT_SQL_ARGS &&\
			jom && nmake clean &&\
			exit\
			" > build.bat

			$COMSPEC //k build.bat
			QT_WAS_BUILT=true
		fi

		export PATH=$PATH:$LIBS_LOCATION/qt/bin

	else
		echo "Qt already installed in $QMAKE_LOCATION"
	fi


	# Install VTK
	cd $LIBS_LOCATION
	if [ ! -d vtk ]; then
		# Download, extract, rename
		download_file http://www.vtk.org/files/release/5.10/$VTK_VERSION.tar.gz ./$VTK_VERSION.tar.gz
		tar -xf $VTK_VERSION.tar.gz
		rm $VTK_VERSION.tar.gz
	# Check for existing installation
	elif [ -f vtk/build/bin/Release/QVTK.lib -a -f vtk/build/bin/Release/vtkRendering.lib ]; then
		if [ $LIB_DEBUG ]; then
			if [ -f vtk/build/bin/Debug/QVTK.lib -a -f vtk/build/bin/Debug/vtkRendering.lib ]; then
				VTK_FOUND=true
			fi
		else
			VTK_FOUND=true
		fi
	fi

	if [ $VTK_FOUND ]; then
		echo "VTK already installed in $LIBS_LOCATION/vtk"
	else
		# Compile
		cd vtk
		mkdir -vp build
		cd build
		cmake .. -DBUILD_TESTING=OFF -DVTK_USE_QT=ON -G "$CMAKE_GENERATOR"
		cmake ..
		cmake --build . --config Release
		cmake --build . --config Release --target QVTK
		cmake --build . --config Debug
		cmake --build . --config Debug --target QVTK
	fi

	# Install shapelib
	cd $LIBS_LOCATION
	if [ ! -d shapelib ]; then
		# Download, extract
		download_file http://download.osgeo.org/shapelib/$SHAPELIB_VERSION.tar.gz ./$SHAPELIB_VERSION.tar.gz
		tar -xf $SHAPELIB_VERSION.tar.gz
		mv $SHAPELIB_VERSION/ shapelib/
		rm -rf $SHAPELIB_VERSION.tar.gz
	elif [ -f shapelib/shapelib.lib ]; then
		SHAPELIB_FOUND=true
	fi

	if [ $SHAPELIB_FOUND ]; then
		echo "Shapelib already installed in $LIBS_LOCATION/shapelib"
	else
		# Compile
		cd shapelib

		echo " \
		\"$WIN_DEVENV_PATH\\..\\..\\VC\\vcvarsall.bat\" $WIN_ARCHITECTURE &&\
		nmake /f makefile.vc &&\
		exit\
		" > build.bat

		$COMSPEC //k build.bat
	fi

	# Install libgeotiff
	cd $LIBS_LOCATION
	if [ ! -d libgeotiff ]; then
		# Download, extract
		download_file http://download.osgeo.org/geotiff/libgeotiff/$LIBGEOTIFF_VERSION.tar.gz ./$LIBGEOTIFF_VERSION.tar.gz
		tar -xf $LIBGEOTIFF_VERSION.tar.gz
		mv $LIBGEOTIFF_VERSION/ libgeotiff/
		rm -rf $LIBGEOTIFF_VERSION.tar.gz
	elif [ -f libgeotiff/geotiff.lib ]; then
		LIBGEOTIFF_FOUND=true
	fi

	if [ $LIBGEOTIFF_FOUND ]; then
		echo "Libgeotiff already installed in $LIBS_LOCATION/libgeotiff"
	else
		# Compile
		cd libgeotiff

		# Download modified makefile
		if [ ! -f makefile_mod.vc ]; then
			download_file http://dl.dropbox.com/u/5581063/makefile_mod.vc ./makefile_mod.vc 14fb13a5bd04ffc298fee7825dc7679f
		fi

		echo " \
		\"$WIN_DEVENV_PATH\\..\\..\\VC\\vcvarsall.bat\" $WIN_ARCHITECTURE &&\
		nmake /f makefile_mod.vc all&&\
		exit\
		" > build.bat

		$COMSPEC //k build.bat
	fi

	# Install Metis
	cd $LIBS_LOCATION
	if [ ! -d metis ]; then
		# Download, extract, rename
		download_file http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/$METIS_VERSION.tar.gz ./$METIS_VERSION.tar.gz
		tar -xf $METIS_VERSION.tar.gz
		mv $METIS_VERSION/ metis/
		rm $METIS_VERSION.tar.gz
	# Check for existing installation
	elif [ -f metis/build/windows/libmetis/Release/metis.lib ]; then
		METIS_FOUND=true
	fi

	if [ $METIS_FOUND ]; then
		echo "Metis already installed in $LIBS_LOCATION/metis"
	else
		# Compile
		cd metis
		$COMSPEC //c "vsgen.bat -G \"$CMAKE_GENERATOR\""
		cd build/windows
		cmake --build . --config Release
		cd ../..
		cp build/windows/libmetis/Release/metis.lib ../lib/metis.lib
		cp include/metis.h ../include/metis.h
	fi

	# Install Boost
	cd $LIBS_LOCATION
	if [ ! -d boost ]; then
		# Download, extract, rename
		download_file http://sourceforge.net/projects/boost/files/boost/$BOOST_VERSION/boost_$BOOST_VERSION_UNDERSCORE.zip/download ./boost_$BOOST_VERSION_UNDERSCORE.zip
		7za x boost_$BOOST_VERSION_UNDERSCORE.zip
		mv boost_$BOOST_VERSION_UNDERSCORE/ boost/
		rm boost_$BOOST_VERSION_UNDERSCORE.zip
	elif [ -d boost/stage/lib ]; then
		BOOST_FOUND=true
	fi

	if [ $BOOST_FOUND ]; then
		echo "Boost is already installed in ..."
	else
		# Compile
		cd boost
		echo " \
			\"$WIN_DEVENV_PATH\\..\\..\\VC\\vcvarsall.bat\" $WIN_ARCHITECTURE &&\
			bootstrap.bat &&\
			bjam.exe &&\
			exit\
			" > build.bat
		$COMSPEC //k build.bat
	fi
fi

cd $SOURCE_LOCATION/scripts/setup