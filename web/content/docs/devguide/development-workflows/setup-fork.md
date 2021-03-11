+++
date = "2021-03-11T11:48"
title = "Set Up your fork"
author = "Lars Bilke"
weight = 1013

aliases = ["/docs/devguide/development-workflows/branching-model"]

[menu]
  [menu.devguide]
    parent = "development-workflows"
+++

<div class='note'>

**Attribution**: The content of this page is largely taken from the [GitHub-blog](https://github.com/blog/2042-git-2-5-including-multiple-worktrees-and-triangular-workflows).
</div>

## Explanation: Forking workflow

Git is very flexible in organizing a distributed development team. We use a so called **Forking workflow**.

The following explanation is taken from an [in-depth article](https://www.atlassian.com/git/tutorials/comparing-workflows#!workflow-forking) on that model:

> Instead of using a single server-side repository to act as the *central* codebase, it gives every developer a server-side repository. This means that each contributor has not one, but two Git repositories: a private local one and a public server-side one.
>
> The main advantage of the Forking Workflow is that contributions can be integrated without the need for everybody to push to a single central repository. Developers push to their own server-side repositories, and only the project maintainer can push to the official repository. This allows the maintainer to accept commits from any developer without giving them write access to the official codebase.
>
> The result is a distributed workflow that provides a flexible way for large, organic teams (including untrusted third-parties) to collaborate securely. This also makes it an ideal workflow for open source projects.
>
> <cite><a href="https://www.atlassian.com/git/tutorials/comparing-workflows#!workflow-forking">www.atlassian.com/git/tutorials/comparing-workflows#!workflow-forking</a> </cite>
>

The workflow is summarized in the following image from the [GitHub blog](https://github.com/blog/2042-git-2-5-including-multiple-worktrees-and-triangular-workflows):
![Git workflow](https://cloud.githubusercontent.com/assets/1319791/8943755/5dcdcae4-354a-11e5-9f82-915914fad4f7.png)

You always **fetch** changes from official repository (called **upstream**), develop on your **local** repository and **push** changes to your server-side repository (called **origin**).

First thing to do when you start working on your local repository is to create a topic branch (based on the current master branch of the official repository) specific to a well defined feature or bugfix you are about to implement. **Never** work on the **master**-branch (it is reserved for the official version)! See also [this tutorial](https://www.atlassian.com/git/tutorials/using-branches) on branching.

Start committing changes in logical chunks. After you are happy with your implementation push your topic branch to your forked repository on GitLab.

Open a [*Merge Request*](https://docs.gitlab.com/ee/user/project/merge_requests/) which will initiate the code review process.

----

## Step: Create a fork

[Create a new fork](https://gitlab.opengeosys.org/ogs/ogs/-/forks/new) from the [official OGS-6 repository](https://gitlab.opengeosys.org/ogs/ogs).

This creates a new fork under your account with the URL `https://gitlab.opengeosys.org/YOUR-USERNAME/ogs`.

## Step: Setup your local clone

You can use the git command line tool to clone the remote repository on GitLab to your PC:

```bash
git clone --filter=blob:limit=100k git@gitlab.opengeosys.org:YOUR-USERNAME/ogs.git
cd ogs
git config remote.pushdefault origin
git config push.default current
```

This creates a new folder `ogs` in your current working directory with the OGS source code. After this step, the remote called `origin` refers to your fork on GitLab. It also sets the default remote for pushes to be `origin` and the default push behavior to `current`. Together this means that if you just type `git push`, the current branch is pushed to the `origin` remote.

<div class='note'>

The `--filter=blob:limit=100k`-parameter instructs git to only fetch files which are smaller than 100 Kbyte. Larger files (e.g. benchmark files, images, PDFs) are fetched on-demand only. This happens automatically and [is a replacement for the previous Git LFS tracked files](https://github.com/ufz/ogs/issues/2961). Requires **git 2.22**!

</div>

Create a second remote called `upstream` that points at the main OGS repository and fetch from it:

```bash
git remote add upstream https://gitlab.opengeosys.org/ogs/ogs.git
git fetch upstream
```

<!-- TODO: rerecord with GitLab -->
<!-- {{< asciinema url="https://asciinema.org/a/249002" speed="3" rows="20" >}} -->

## Optional: Enable git commit hooks

[Git hooks](https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks) help to check for several issues before doing commits or pushes and it is highly recommended to enable these checks.

Install [pre-commit](https://pre-commit.com/) (a git hook manager) via Pythons `pip`:

```bash
pip3 install --user pre-commit
```

This installed `pre-commit` to `.local/bin` in your home directory or to `C:\Users\[username]\AppData\Roaming\Python\Python37\Scripts` on Windows. Make sure to have this directory in your `PATH`!

Enable the hooks in the source code with:

```bash
cd ogs
pre-commit install
```

## Optional: Working on a new feature

You only have to follow the above steps once. From then on, whenever you want to work on a new feature, you can more easily interact with the remote repositories.

Make sure that your local repository is up-to-date with the upstream repository:

```bash
git fetch upstream
```

Create a branch `feature-name` off of upstream `master`-branch to work on a new feature, and check out the branch:

```bash
git checkout -b feature-name upstream/master
```

This automatically sets up your local `new-feature`-branch to track the upstream `master`-branch. This means that if more commits are added to `master` upstream, you can merge those commits into your `feature`-branch by typing

```bash
git pull
```

or rebase your branch on top of the new master by typing

```bash
git pull --rebase
```

Now after you implemented the feature and committed your work you can push the new commits to the `feature-name`-branch on your GitLab fork:

```bash
git push
```

If your work is done submit a [merge request](https://gitlab.opengeosys.org/ogs/ogs/-/merge_requests/new).

This workflow is summarized with this picture:
![Workflow visualization](https://cloud.githubusercontent.com/assets/1319791/8943755/5dcdcae4-354a-11e5-9f82-915914fad4f7.png)
