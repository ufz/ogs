DIR=$( cd $(dirname "${BASH_SOURCE[0]}") ; pwd -P )
source $DIR/cli.sh
export CC=`which mpicc`
export CXX=`which mpicxx`
