if [ -n "$ZSH_VERSION" ]; then
    DIR=$( cd $(dirname "${(%):-%x}") ; pwd -P )
else
    DIR=$( cd $(dirname "${BASH_SOURCE[0]}") ; pwd -P )
fi

source $DIR/cli.sh
source $DIR/user.sh

module load HDF5/1.14.0
module load Python/3.10.8

# PETSc and VTK build by ogs