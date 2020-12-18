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
- Create new web release page with generated artifact names and changelog (convert MR ids to urls: replace `!([0-9][0-9][0-9][0-9])` with `[!$1](https://gitlab.opengeosys.org/ogs/ogs/-/merge_requests/$1)` and `#([0-9][0-9][0-9][0-9])` with `[#$1](https://gitlab.opengeosys.org/ogs/ogs/-/issues/$1)`)
- Create a tag and push
- A new release is automatically created on GitLab
  - Fill in the release notes from the Wiki
- Copy release binaries and container images from CI job to Azure OGS storage to a subdirectory containing the tag name at <https://ogsstorage.blob.core.windows.net/binaries/ogs6>
- Create a release on GitHub mirror (`ufz/ogs`)
- Check if a [Zenodo release](https://zenodo.org/account/settings/github/repository/ufz/ogs#) is automatically issued
- Issue a scan on [Software Heritage Archive](https://archive.softwareheritage.org/browse/origin/directory/?origin_url=https://gitlab.opengeosys.org/ogs/ogs.git)
- Create a discourse announcement post
