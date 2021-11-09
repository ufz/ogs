module load CMake/3.18.0
module load GCC/10.3.0
module load ParaStationMPI/5.4.10-1
module load PETSc/3.14

export CC=mpicc
export CXX=mpic++

echo "Example config: cmake ../ogs -DCMAKE_BUILD_TYPE=Release -DOGS_USE_PETSC=ON"
echo "Example run: srun bin/ogs -n 3 -t 2 -A ogs6hpc3 bin/ogs ../ogs/Tests/Data/EllipticPETSc/cube_1e3_neumann.prj"
