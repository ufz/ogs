+++
date = "2018-02-23T15:28:13+01:00"
title = "Get the source code"
author = "Lars Bilke"
weight = 3

[menu]
  [menu.devguide]
    parent = "getting-started"
+++

<div class='note'>

<i class="far fa-exclamation-triangle"></i> This page describes how to get the source code as a simple download or git clone. The information on the forking workflow with `git` has been moved to the [Development workflows]({{< ref "setup-fork.md" >}})-section.

</div>

## Clone the source code repository with Git

First you need to get the clone URL:

![Getting the git clone URL from the GitLab web interface](git-url.png)

Then clone the repository with `git`:

```bash
git clone --filter=blob:limit=100k https://gitlab.opengeosys.org/ogs/ogs.git
```

<div class='note'>

The `--filter=blob:limit=100k`-parameter instructs git to only fetch files which are smaller than 100 Kilobyte. Larger files (e.g. benchmark files, images, PDFs) are fetched on-demand only. This happens automatically and [is a replacement for the previous Git LFS tracked files](https://gitlab.opengeosys.org/ogs/ogs/-/issues/2961). Requires at least **git 2.22**!

</div>
