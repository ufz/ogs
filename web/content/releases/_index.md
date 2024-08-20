+++
Title = "Releases"

[menu.main]
name = "Releases"
weight = 1

aliases = ["/downloads", "/download"] # Redirect for Hydrology II tutorial

[[head_downloads]]
name = "Latest pip-installable package"
url = "/docs/userguide/basics/introduction#installation"
note = "All platforms: `pip install --pre --index-url https://gitlab.opengeosys.org/api/v4/projects/120/packages/pypi/simple ogs`"
icon = "far fa-arrow-right"

# url encode job name, https://www.w3schools.com/tags/ref_urlencode.ASP
[[head_downloads]]
name = "Latest Windows CLI"
url = "https://gitlab.opengeosys.org/ogs/ogs/-/jobs/artifacts/master/browse/build?job=build+win"
note = "Download and unpack .zip-file"
icon = "fab fa-windows"

[[head_downloads]]
name = "Latest Windows Data Explorer"
url = "https://gitlab.opengeosys.org/ogs/ogs/-/jobs/artifacts/master/browse/build?job=build+gui+win"
note = "Download and unpack .zip-file"
icon = "fab fa-windows"

[[head_downloads]]
name = "Latest container CLI"
url = "https://minio.ufz.de/ogs/public/container/ogs/master/ogs-serial.squashfs"
note = "Download and run .squashfs-file with [Singularity / Apptainer](/docs/userguide/basics/container/)"
icon = "far fa-container-storage"

[[head_downloads]]
name = "Latest container CLI with PETSc-support"
url = "https://minio.ufz.de/ogs/public/container/ogs/master/ogs-petsc.squashfs"
note = "Download and run .squashfs-file with [Singularity / Apptainer](/docs/userguide/basics/container/)"
icon = "far fa-container-storage"

[[head_downloads]]
name = "Latest container CLI with Intel MKL-support"
url = "https://minio.ufz.de/ogs/public/container/ogs/master/ogs-mkl.squashfs"
note = "Download and run .squashfs-file with [Singularity / Apptainer](/docs/userguide/basics/container/)"
icon = "far fa-container-storage"

[[head_downloads]]
name = "Latest container CLI with PETSc- and Intel MKL-support"
url = "https://minio.ufz.de/ogs/public/container/ogs/master/ogs-petsc-mkl.squashfs"
note = "Download and run .squashfs-file with [Singularity / Apptainer](/docs/userguide/basics/container/)"
icon = "far fa-container-storage"

[[head_downloads]]
name = "Instructions on downloading latest benchmark input files"
url = "/docs/userguide/basics/introduction/#download-benchmarks"
icon = "far fa-arrow-right"
+++
