+++
date = "2020-05-29T11:00:13+01:00"
title = "Migrate to GitLab"
author = "Lars Bilke"
weight = 1043

[menu]
  [menu.devguide]
    parent = "advanced"
+++

To migrate your repository run the following commands in your local repository (`cd path/to/ogs`):

## Disable Git LFS

We used Git LFS to store large files but [switched back to plain git](https://github.com/ufz/ogs/issues/2961).

<div class='note'>

### <i class="far fa-info-circle"></i> IMPORTANT: Normalize LFS files in your existing branches

If you have branches from pre-GitLab times with e.g. newly created LFS files (benchmark files, images, ...) you have to convert them back to plain git:

```bash
git add --renormalize .
git commit -m "Converted LFS files to plain git."
```

Now you have to squash this conversion commit into your original commit which added the files as Git LFS files. In result your branch history should and **must not** have any Git LFS files! Otherwise GitLab will reject the push!
</div>

When you are done migrating your branches you need to disable Git LFS in your local repo:

```bash
git lfs uninstall --local
```

## Create a new fork on GitLab

[Create a new fork](https://gitlab.opengeosys.org/ogs/ogs/-/forks/new) from the [official OGS-6 repository](https://gitlab.opengeosys.org/ogs/ogs).

This creates a new fork under your account with the URL `https://gitlab.opengeosys.org/YOUR-USERNAME/ogs`.

## Migrate your local repos to point to GitLab

You have to modify your git remotes to point to the new GitLab repos. Assuming the former official git repo remote is called `upstream`:

```bash
git remote set-url upstream https://gitlab.opengeosys.org/ogs/ogs.git
```

Assuming your personal forks remote is called `origin`:

```bash
git remote set-url origin git@gitlab.opengeosys.org:YOUR-USERNAME/ogs.git
```

----

Or you can clone a fresh repo by following the steps in [Get the source code]({{< ref "get-the-source-code.md" >}}).
