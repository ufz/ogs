# previous version Stages/2023 can be used as long as it is supported by JUWELS
# the module versions can be found in git history
module load Stages/2024
module load GCC
module load ParaStationMPI
module load CMake
module load Ninja
module load Python
module load Eigen
module load Boost
module load git
module load PETSc
module load netCDF-C++4
# HDF5 version 1.14.3 fixes a bug that is contained in previous versions
# a module with this version isn't in Stages/2024
# it will be provided by ogs cpm
# module load HDF5/1.14.3
module load imkl # provides lapack - use either imkl or ScaLAPACK module
# module load ScaLAPACK

export CC=mpicc
export CXX=mpic++

# CMake may spawn as many processes as available. This may overload the filesystem on JUWELS.
# It is therefore necessary to limit the number of processes with CMAKE_BUILD_PARALLEL_LEVEL.
export CMAKE_BUILD_PARALLEL_LEVEL=6

echo "Example config: cmake --preset release-petsc"
echo "Example run: srun -n 3 -t 2 -A ogs6hpc4 bin/ogs ../ogs/Tests/Data/EllipticPETSc/cube_1e3_neumann.prj"
