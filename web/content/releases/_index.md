+++
Title = "Releases"

[menu.main]
name = "Releases"
weight = 1

aliases = ["/downloads", "/download"] # Redirect for Hydrology II tutorial

[[head_downloads]]
name = "Latest pip-installable package"
url = "/docs/userguide/basics/introduction#get-current-development-version-with-pip"
note = "All platforms: `pip install --pre --index-url https://gitlab.opengeosys.org/api/v4/projects/120/packages/pypi/simple ogs`"
icon = "far fa-arrow-right"

# url encode job name, https://www.w3schools.com/tags/ref_urlencode.ASP
[[head_downloads]]
name = "Latest Windows CLI with Utilities"
url = "https://gitlab.opengeosys.org/ogs/ogs/-/jobs/artifacts/master/browse/build/release?job=build+win"
note = "Download and unpack .zip-file"
icon = "fab fa-windows"

[[head_downloads]]
name = "Latest Windows Data Explorer with Utilities"
url = "https://gitlab.opengeosys.org/ogs/ogs/-/jobs/artifacts/master/browse/build/release-gui?job=build+gui+win"
note = "Download and unpack .zip-file"
icon = "fab fa-windows"

[[head_downloads]]
name = "Latest Singularity container CLI with Utilities"
url = "https://gitlab.opengeosys.org/ogs/ogs/-/jobs/artifacts/master/browse/ThirdParty/container-maker/_out/images?job=container"
note = "Download and run .sif-file with [Singularity](/docs/userguide/basics/container/)"
icon = "far fa-container-storage"

[[head_downloads]]
name = "Instructions on downloading latest benchmark input files"
url = "/docs/userguide/basics/introduction/#download-benchmarks"
icon = "far fa-arrow-right"
+++
