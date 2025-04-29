if [ -n "$ZSH_VERSION" ]; then
    DIR=$( cd $(dirname "${(%):-%x}") ; pwd -P )
else
    DIR=$( cd $(dirname "${BASH_SOURCE[0]}") ; pwd -P )
fi

source $DIR/cli.sh
export CC=mpicc
export CXX=mpic++

# Fixes:
# Error obtaining unique transport key from PMIX
# (OMPI_MCA_orte_precondition_transports not present in the environment).
# https://users.open-mpi.narkive.com/p7v4NTFg/ompi-error-when-calling-mpi-init
eval "export `mpirun env | grep OMPI_MCA_orte_precondition_transports | head -n 1`"
