if [ ! -f ~/.easybuild-yes ]; then
    echo "ERROR: Easybuild modules not enabled but required!\n"
    echo "Run 'touch ~/.easybuild-yes' and re-login to enable."
    echo "For more details see:\n  https://www.opengeosys.org/docs/devguide/advanced/working-on-eve"
    return 1
fi

module use /global/apps/modulefiles

module load foss/2019b
module load CMake/3.15.3
module load ninja
module load git/2.23.0
module load git-lfs/2.7.1

# Libraries
module load boost/1.67.0-1
module load eigen/3.3.4-1-cmake
module load vtk/8.2.0/foss2019b/serial
module load HDF5/1.10.5-nompi

# Tools
module load ccache/3.3.3
