export CMAKE_BUILD_PARALLEL_LEVEL="${CMAKE_BUILD_PARALLEL_LEVEL:-8}"
export CTEST_PARALLEL_LEVEL="${CTEST_PARALLEL_LEVEL:-4}"

module use /global/apps/modulefiles

module load foss/2024a
module load CMake/3.29.3
module load Ninja/1.12.1

# Tools
module load ccache/3.3.3
module load git-lfs
module load Python/3.12.3

# Python dependencies
python -m venv .venv
source .venv/bin/activate
pip install "numpy<2"
