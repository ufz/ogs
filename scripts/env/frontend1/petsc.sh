if [ -n "$ZSH_VERSION" ]; then
    DIR=$( cd $(dirname "${(%):-%x}") ; pwd -P )
else
    DIR=$( cd $(dirname "${BASH_SOURCE[0]}") ; pwd -P )
fi

source $DIR/mpi.sh
module load petsc/3.7.6_maint_gcc6.2.0_openmpi_gcc_1.8.8-1

echo -e "Note: If you want to run a simulation on the cluster be aware of the "\
     "mixed CPU architecture. There are Sandy-Bridge-based nodes (orte-28, "\
     "frontend2, envinf1) as well as Skylake-based nodes (orte-40, frontend1)"\
     ".\nConsider setting CMake-option OGS_CPU_ARCHITECTURE to:\n\n"\
     "  -DOGS_CPU_ARCHITECTURE=sandybridge\n"\
     "   -- OR --\n"\
     "  -DOGS_CPU_ARCHITECTURE=skylake-avx512"
