+++
date = "2018-02-26T11:00:13+01:00"
title = "Publish a release"
author = "Lars Bilke"
weight = 1080

[menu]
  [menu.devguide]
    parent = "procedures"
+++

## Publication procedure

- Update merge request template (settings / merge_requests) to point to a new changelog wiki page
- Run `python scripts/python/do-release.py`
- Update `CITATION.cff`, create a commit, tag and push (see script output)
- A new release is automatically created on GitLab
  - Fill in the release notes from the Wiki
  - Convert MR ids to URLs: replace `!([0-9][0-9][0-9][0-9])` with `[!$1](https://gitlab.opengeosys.org/ogs/ogs/-/merge_requests/$1)` and `#([0-9][0-9][0-9][0-9])` with `[#$1](https://gitlab.opengeosys.org/ogs/ogs/-/issues/$1)`
- Copy release binaries and container images from CI job to Azure OGS storage to a subdirectory containing the tag name at <https://ogsstorage.blob.core.windows.net/binaries/ogs6>
- Create a release on GitHub mirror (`ufz/ogs`)
- Check if a [Zenodo release](https://zenodo.org/account/settings/github/repository/ufz/ogs#) is automatically issued
- Issue a scan on Software Heritage: `curl -X POST https://archive.softwareheritage.org/api/1/origin/save/git/url/https://gitlab.opengeosys.org/ogs/ogs.git/`
- Update `CITATION.cff` and `web/content/publications/_index.md` with new Zenodo DOI
- Update `CITATION.cff` author list (`git shortlog -sne 6.4.3...6.4.2`) and corresponding bibtex-entry in publications web page, TODO: download BibTex entry and copy to web page, e.g. https://zenodo.org/records/2600042/export/bibtex
- Update Zenodo entry with correct authors
- Create a discourse announcement post
  - DOI badge
  - Software heritage badge
  - Link to release page on opengeosys.org
  - Link to GitLab release
  - Link to Zenodo release
  - Link to archive on Software Heritage
- Set milestone on expired issues and merge requests
