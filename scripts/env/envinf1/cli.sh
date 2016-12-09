module () { eval `/usr/local/modules/3.2.10-1/Modules/3.2.10/bin/modulecmd sh $*`; }
export MODULEPATH=$MODULEPATH:/global/apps/modulefiles

module load cmake/3.1.3-1
module load gcc/4.8.1-3

# Libraries
module load boost/1.62.0-1
module load eigen/3.2.8-1-cmake
module load vtk/6.3.0_openmpi-1.8.4-noqt-1

# Tools
module load numdiff/5.8.1-1
module load coreutils/8.21-1
