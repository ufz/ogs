if [ -n "$ZSH_VERSION" ]; then
    DIR=$( cd $(dirname "${(%):-%x}") ; pwd -P )
else
    DIR=$( cd $(dirname "${BASH_SOURCE[0]}") ; pwd -P )
fi

source $DIR/mpi.sh

# TODO build these for foss2022b or remove it:
#module load petsc/3.18.5_foss2020b
#module load vtk/9.2.5_mpi_foss2020b
#module load hdf5/1.14.0_mpi_foss_2020b

export OGS_PETSC_CONFIG_OPTIONS="--download-fc;--download-mumps;--download-hypre;--download-scalapack"

echo -e "Note: If you want to run a simulation on the cluster be aware of the"\
     "mixed CPU architecture. There are Sandy-Bridge-based nodes (orte-28,"\
     "frontend2) as well as Skylake-based nodes (orte-40, frontend1)"\
     ".\nConsider setting CMake-option OGS_CPU_ARCHITECTURE to:\n\n"\
     "  -DOGS_CPU_ARCHITECTURE=sandybridge\n"\
     "   -- OR --\n"\
     "  -DOGS_CPU_ARCHITECTURE=skylake-avx512"
