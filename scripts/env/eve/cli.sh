if [ ! -f ~/.easybuild-yes ]; then
    echo "ERROR: Easybuild modules not enabled but required!\n"
    echo "Run 'touch ~/.easybuild-yes' and re-login to enable."
    echo "For more details see:\n  https://www.opengeosys.org/docs/devguide/advanced/working-on-eve"
    return 1
fi

module use /global/apps/modulefiles

module load foss/2018b
module load cmake
module load ninja
module load git-lfs

# Libraries
module load Boost/1.67.0
module load eigen/3.3.4-1-cmake
module load vtk/8.2.0/foss2018b/serial

# Tools
module load ccache/3.3.3
