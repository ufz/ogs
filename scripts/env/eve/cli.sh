export CMAKE_BUILD_PARALLEL_LEVEL="${CMAKE_BUILD_PARALLEL_LEVEL:-8}"
export CTEST_PARALLEL_LEVEL="${CTEST_PARALLEL_LEVEL:-4}"

module use /global/apps/modulefiles

module load foss/2022b
module load cmake/3.22.4-1
module load Ninja/1.11.1

# Tools
module load ccache/3.3.3
module load git-lfs
module load Python/3.10.8

# Python dependencies
virtualenv .venv
source .venv/bin/activate
pip install numpy
