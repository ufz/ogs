## ogs-data ##

WARNING: This repository is a mirror of test data from the
[OGS source code](https://github.com/ufz/ogs/tree/master/Tests/Data).

It is provided for convenience and can be used with a precompiled `ogs`
executable or a OGS Singularity container.

### Usage

Run all benchmarks with [bats](https://github.com/bats-core/bats-core) in a
singularity container:

```bash
OGS_VERSION=6.2.2
curl -L -O https://github.com/ufz/ogs/releases/download/$OGS_VERSION/ogs-$OGS_VERSION-serial.sif
export SIF=ogs-$OGS_VERSION-serial.sif
git clone --depth=1 --branch $OGS_VERSION https://github.com/ufz/ogs-data.git
export SRC=$PWD/ogs-data
export OUT=$PWD/out
# Note: the .bats file is included since OGS > 6.2.2, you can try to copy it
# into older source trees
bats ogs-data/benchmarks.bats
```
