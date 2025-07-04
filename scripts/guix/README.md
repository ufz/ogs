
## Runtime

```bash
# builds ogs serial config and starts isolated shell (like in a container)
guix time-machine -C scripts/guix/channels.scm -- shell -C -m scripts/guix/manifest.scm
# ogs petsc config
guix time-machine -C scripts/guix/channels.scm -- shell -C -m scripts/guix/manifest-petsc.scm
# To create an archivable Apptainer container:
guix time-machine -C scripts/guix/channels.scm -- pack -m scripts/guix/manifest.scm \
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
