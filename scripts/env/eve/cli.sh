if [ ! -f ~/.easybuild-yes ]; then
    echo "ERROR: Easybuild modules not enabled but required!\n"
    echo "Run 'touch ~/.easybuild-yes' and re-login to enable."
    echo "For more details see:\n  https://www.opengeosys.org/docs/devguide/advanced/working-on-eve"
    return 1
fi

module use /global/apps/modulefiles

module load cmake
module load foss/2018b
module load ninja/1.9.0
module load git/2.20.1

# Libraries
module load boost/1.62.0-1
module load eigen/3.3.4-1-cmake
module load vtk/8.2.0/foss2018b/serial

# Tools
module load coreutils/8.21-1
module load ccache/3.3.3
