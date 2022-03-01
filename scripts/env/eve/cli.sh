if [ ! -f ~/.easybuild-yes ]; then
    echo "ERROR: Easybuild modules not enabled but required!\n"
    echo "Run 'touch ~/.easybuild-yes' and re-login to enable."
    echo "For more details see:\n  https://www.opengeosys.org/docs/devguide/advanced/working-on-eve"
    return 1
fi

module use /global/apps/modulefiles

module load foss/2020b
module load cmake/3.21.1-1
module load Ninja/1.10.1

# Tools
module load ccache/3.3.3
module load git-lfs
module load Python/3.8.6
