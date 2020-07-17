+++
date = "2018-02-26T11:00:13+01:00"
title = "GitLab CI"
author = "Lars Bilke"
weight = 1022

[menu]
  [menu.devguide]
    parent = "testing"
+++

## Introduction

[GitLab CI](https://docs.gitlab.com/ee/ci/) is a powerful [Continuous Integration](../../development-workflows/continuous-integration) system integrated into GitLab.

The tasks of the CI system are configured in [scripts inside the OGS source code](https://gitlab.opengeosys.org/ogs/ogs/-/tree/master/scripts/ci). The entrypoint is defined in [.gitlab-ci.yml](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/.gitlab-ci.yml). Scripting and versioning the configuration together with the source code is very powerful, e.g. if you introduce a new OGS CMake configuration in a merge request even the change of the CI jobs configuration or jobs environment (Docker container definition) can be part of the merge request.

## GitLab Pipeline

A CI run consists of a [pipeline](https://docs.gitlab.com/ee/ci/pipelines/) which contains [stages](https://docs.gitlab.com/ee/ci/yaml/#stages) which in turn contain jobs. A job runs a set of instructions (e.g. checking out the source code, building the code, testing the code) on a [runner](https://docs.gitlab.com/runner/).

Each pipeline run is visualized as follows:
![GitLab pipeline visualization](../GitLab-Pipeline.png)

Jobs are belong to a stage and each job will get a status (success, warnings, failure). Some jobs are optional (see the gear icon) and can be manually triggered by pressing the play button.

## Automatic testing

The master-branch of the the main repository as well as all merge requests on that repo are automatically tested. See [the pipelines page](https://gitlab.opengeosys.org/ogs/ogs/pipelines).

### Skip automatic testing

If you want to skip a pipeline run for a push add the `-o ci.skip` git push option. Example:

```bash
git push -o ci.skip
```

Or add add `[ci skip]` to the commit message to skip the pipeline for this commit. Example:

```bash
git commit -m "Added feature X [ci skip]"
```

### Automatic testing for own repository

TODO
