
## Runtime

```bash
# builds ogs serial config and starts isolated shell (like in a container)
guix time-machine -C scripts/guix/channels.scm -- shell -C -m scripts/guix/manifest-serial.scm
# ogs petsc config
guix time-machine -C scripts/guix/channels.scm -- shell -C -m scripts/guix/manifest-petsc.scm
# To create an archivable Apptainer container:
guix time-machine -C scripts/guix/channels.scm -- pack -m scripts/guix/manifest-serial.scm \
  -RR -f squashfs
```

## Development

```bash
# Start development shell:
guix time-machine -C scripts/guix/channels.scm -- shell -C -m scripts/guix/manifest-dev.scm

# To create an archivable Apptainer container:
guix time-machine -C scripts/guix/channels.scm -- pack -f squashfs -m scripts/guix/manifest-dev.scm
...
apptainer shell  /gnu/store/...-bash-glibc-locales-nss-certs-coreutils-squashfs-pack.gz.squashfs
# Now in the container
cmake -G Ninja `realpath /gnu/store/*-ogs-source-1.0` -B build --preset release # or debug
cd build
ninja
ninja tests
ctest -LE large -j 16
```

## Run ctests with Guix OGS executables

To run ctests in Guix, set `OGS_TEST_USE_PATH` to `ON`. CTest then resolves
project executables from the Guix shell's `PATH`.

`guix shell` caches the profile of a manifest file and reuses it until the
manifest file changes, so source changes would not rebuild OGS;
`--rebuild-cache` forces the rebuild.

Usage:

```bash
# Run from the repository root on the host.
cmake --preset release -DOGS_TEST_USE_PATH=ON \
      -DOGS_EXCLUDE_CTESTS='xmlpatch;partmesh' -B ../build/guix-test

# Run CTest against the same build directory, within the Guix environment.
guix time-machine -C scripts/guix/channels.scm -- \
  shell --rebuild-cache -m scripts/guix/manifest-serial.scm cmake -- \
  ctest --test-dir ../build/guix-test -LE large --output-on-failure -j 8
```

The Guix OGS build does not install the repository's `apply_xmlpatch` test
helper, and its `partmesh` tests are excluded for this workflow.
