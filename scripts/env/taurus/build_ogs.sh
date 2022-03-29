# Alternative: OGS_DIR=$HOME/o/
OGS_DIR=$( realpath ../../../.. )
BIN_DIR=rp
ENV=taurus
source $OGS_DIR/s/scripts/env/$ENV/petsc.sh

echo "Folder with all PETSC binaries for selected Partition: $PETSC_DIR"

BUILD_DIR=$OGS_DIR/build/$PARTITION_NAME/$BIN_DIR/

mkdir -p $OGS_DIR/build/cpm_cache
mkdir -p $BUILD_DIR
cd $BUILD_DIR

export CPM_SOURCE_DIR=$OGS_DIR/build/cpm_cache

CC=`which mpicc` CXX=`which mpic++` cmake $OGS_DIR/s -G Ninja -DCMAKE_BUILD_TYPE=Release -DOGS_EIGEN_DYNAMIC_SHAPE_MATRICES=Off -DOGS_USE_PCH=Off -DOGS_USE_PETSC=On -DOGS_BUILD_PROCESSES="HT" -DPETSC_DIR=$PETSC_DIR -DPETSC_ARCH=arch-linux-c-opt

ninja -j 8