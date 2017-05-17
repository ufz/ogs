DIR=$( cd $(dirname "${BASH_SOURCE[0]}") ; pwd -P )
source $DIR/mpi.sh
module load petsc-bilke/3.7.6_gcc-6.2.0_openmpi-1.8.8
