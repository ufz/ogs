module () { eval `/usr/local/modules/3.2.10-1/Modules/3.2.10/bin/modulecmd sh $*`; }
export MODULEPATH=$MODULEPATH:/global/apps/modulefiles

module load python/2
module load cmake/3.6.2-1
module load gcc/6.2.0-1

# Tools
module load coreutils/8.21-1
module load ccache/3.3.3

source /global/apps/ogs/virtualenv/conan/bin/activate
