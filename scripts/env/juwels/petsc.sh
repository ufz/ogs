module load Stages/2023
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load CMake/3.23.1
module load Ninja/1.10.2
module load Python
module load Eigen
module load Boost
module load git
module load HDF5/1.12.2
module load ScaLAPACK/2.2.0-fb

export CC=mpicc
export CXX=mpic++

echo "Example config: cmake --preset release-petsc"
echo "Example run: srun -n 3 -t 2 -A ogs6hpc4 bin/ogs ../ogs/Tests/Data/EllipticPETSc/cube_1e3_neumann.prj"
