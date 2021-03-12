if [ ! -f ~/.easybuild-yes ]; then
    echo "ERROR: Easybuild modules not enabled but required!\n"
    echo "Run 'touch ~/.easybuild-yes' and re-login to enable."
    echo "For more details see:\n  https://www.opengeosys.org/docs/devguide/advanced/working-on-eve"
    return 1
fi

module use /global/apps/modulefiles

module load foss/2019b
module load cmake/3.19.4-1
module load ninja
module load git/2.23.0

# Tools
module load ccache/3.3.3
