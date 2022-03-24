if [ -n "$ZSH_VERSION" ]; then
    DIR=$( cd $(dirname "${(%):-%x}") ; pwd -P )

else
    DIR=$( cd $(dirname "${BASH_SOURCE[0]}") ; pwd -P )
fi

source $DIR/petsc.sh

PETSC_INSTALL=$ALL_PETSC_DIR/$PARTITION_NAME/
mkdir -p $PETSC_INSTALL

cd $ALL_PETSC_DIR
echo "Clone or checkout PETSC (Version: $PETSC_VERSION) to $PETSC_DIR."
git clone --depth 1 --branch $PETSC_VERSION --single-branch https://gitlab.com/petsc/petsc.git $PETSC_DIR || (cd $PETSC_DIR ; git pull)
mkdir -p petsc_$PETSC_VERSION
cd $PETSC_DIR
echo "Build PETSC in $PETSC_DIR started."
./configure --with-cc=gcc --with-cxx=g++ --with-fc=0 --prefix=PETSC_INSTALL --with-debugging=0
make
make PETSC_DIR=$PETSC_DIR PETSC_ARCH=arch-linux-c-opt install
echo "Build PETSC in $PETSC_DIR done."
