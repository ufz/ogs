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
- Copy release binaries and container images from CI job to OGS S3 storage to a subdirectory containing the tag name at <https://vip.s3.ufz.de/ogs/public/binaries/ogs6>
  - Create folder `mc mb ogs/ogs/public/binaries/ogs6/TAG/`
  - Copy `mc cp ogs/ogs/public/container/ogs/master/* ogs/ogs/public/binaries/ogs6/TAG/`
- Create a release on GitHub mirror (`ufz/ogs`)
- Check if a [Zenodo release](https://zenodo.org/account/settings/github/repository/ufz/ogs#) is automatically issued
- Update Zenodo entry with correct authors (obtained via `git shortlog -sne [new_version]...[previous_version]`)
- Run `python scripts/python/post-release.py` and commit and create a discourse announcement post
- Set milestone on expired issues and merge requests
