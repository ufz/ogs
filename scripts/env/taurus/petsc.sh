if [ -n "$ZSH_VERSION" ]; then
    DIR=$( cd $(dirname "${(%):-%x}") ; pwd -P )
else
    DIR=$( cd $(dirname "${BASH_SOURCE[0]}") ; pwd -P )
fi

source $DIR/cli.sh
source $DIR/user.sh

module load HDF5/1.12.1
module load Python/3.9.6
module load Boost/1.77

# VTK will be built by CPM, because no VTK build with GCCCore/11.2.0 in foss/2021b available

# PETSC needs custom build!
PETSC_VERSION=v3.16.4
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ALL_PETSC_DIR/$PARTITION_NAME/lib
PETSC_DIR=$ALL_PETSC_DIR/$PARTITION_NAME/petsc_$PETSC_VERSION