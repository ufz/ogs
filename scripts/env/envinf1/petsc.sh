DIR=$( cd $(dirname "${BASH_SOURCE[0]}") ; pwd -P )
source $DIR/mpi.sh
module load petsc/3.7.6_maint_gcc6.2.0_openmpi_gcc_1.8.8-1
