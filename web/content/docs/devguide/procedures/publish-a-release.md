+++
date = "2018-02-26T11:00:13+01:00"
title = "Publish a release"
author = "Lars Bilke"
weight = 1051

[menu]
  [menu.devguide]
    parent = "procedures"
+++

## Publication procedure

- Update [merge request template](https://gitlab.opengeosys.org/ogs/ogs/edit) to point to a new changelog wiki page
- Update `CHANGELOG.md` to point to new GitLab release
- Create new web release page with generated artifacts
- Create a tag
- Create a [new release on GitLab](https://gitlab.opengeosys.org/ogs/ogs/-/tags/new)
  - Fill in the message, e.g. "OpenGeoSys Release 6.0.8"
  - Fill in the release notes from the Wiki
- Copy release binaries and container images from CI job to Azure OGS storage to a subdirectory containing the tag name at <https://ogsstorage.blob.core.windows.net/binaries/ogs6>
- Check if a [Zenodo release](https://zenodo.org/account/settings/github/repository/ufz/ogs#) is automatically issued
- Issue a scan on [Software Heritage Archive](https://archive.softwareheritage.org/browse/origin/directory/?origin_url=https://gitlab.opengeosys.org/ogs/ogs.git)
- Create a discourse announcement post
