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
- Create new netlify site (in an empty directory)
  <!-- vale off -->
  - `netlify init`
  - `# [ENTER]`
  - `# ogs-doxygen-[TAG (- separated instead of .)]`
  <!-- vale on -->
- Update `CITATION.cff`, create a commit, tag and push (see script output)
- A new release is automatically created on GitLab
  - Fill in the release notes from the Wiki
  - Convert MR ids to URLs: replace `!([0-9][0-9][0-9][0-9])` with `[!$1](https://gitlab.opengeosys.org/ogs/ogs/-/merge_requests/$1)` and `#([0-9][0-9][0-9][0-9])` with `[#$1](https://gitlab.opengeosys.org/ogs/ogs/-/issues/$1)`
- Copy release binaries and container images from CI job to Azure OGS storage to a subdirectory containing the tag name at <https://ogsstorage.blob.core.windows.net/binaries/ogs6>
- Create a release on GitHub mirror (`ufz/ogs`)
- Check if a [Zenodo release](https://zenodo.org/account/settings/github/repository/ufz/ogs#) is automatically issued
- Run `python scripts/python/post-release.py` and commit and create a discourse announcement post
- Update Zenodo entry with correct authors (obtained via `git shortlog -sne [new_version]...[previous_version]`)
- Set milestone on expired issues and merge requests
