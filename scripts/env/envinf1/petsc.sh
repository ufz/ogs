DIR=$( cd $(dirname "${BASH_SOURCE[0]}") ; pwd -P )
source $DIR/mpi.sh
module load petsc/3.7.2_maint_petsc_maint_3.7.2_gcc_openmpi_1.8.4-2
