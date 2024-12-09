+++
date = "2021-03-11T11:48:00"
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

<!-- vale off -->

> Instead of using a single server-side repository to act as the *central* codebase, it gives every developer a server-side repository. This means that each contributor has not one, but two Git repositories: a private local one and a public server-side one.
>
> The main advantage of the Forking Workflow is that contributions can be integrated without the need for everybody to push to a single central repository. Developers push to their own server-side repositories, and only the project maintainer can push to the official repository. This allows the maintainer to accept commits from any developer without giving them write access to the official codebase.
>
> The result is a distributed workflow that provides a flexible way for large, organic teams (including not trusted third-parties) to collaborate securely. This also makes it an ideal workflow for open source projects.
>
> <cite><a href="https://www.atlassian.com/git/tutorials/comparing-workflows#!workflow-forking">www.atlassian.com/git/tutorials/comparing-workflows#!workflow-forking</a> </cite>
>

<!-- vale on -->

The workflow is summarized in the following image from the [GitHub blog](https://github.com/blog/2042-git-2-5-including-multiple-worktrees-and-triangular-workflows):
![Git workflow](https://cloud.githubusercontent.com/assets/1319791/8943755/5dcdcae4-354a-11e5-9f82-915914fad4f7.png)

<div class='note'>

### Important naming conventions

You always **fetch** changes from the official repository (called `upstream`), develop on your **local** repository and **push** changes to your personal server-side repository (called `origin`).

</div>

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
```

This creates a new folder `ogs` in your current working directory with the OGS source code. After this step, the remote called `origin` refers to your fork on GitLab.

<div class='note'>

The `--filter=blob:limit=100k`-parameter instructs git to only fetch files which are smaller than 100 Kilobyte. Larger files (e.g. benchmark files, images, PDFs) are fetched on-demand only. This happens automatically and [is a replacement for the previous Git LFS tracked files](https://gitlab.opengeosys.org/ogs/ogs/-/issues/2961). Requires at least **git 2.22**!

</div>

### Local clone git settings

After that initially set some useful git settings for your local repo:

```bash
cd ogs
```

The following sets the default remote for pushes to be `origin` and the default push behavior to `current`. Together this means that if you just type `git push`, the current branch is pushed to the `origin` remote.

```bash
git config remote.pushdefault origin
git config push.default current
```

To streamline the updating workflow the `rebase`-command is configured to handle local modifications automatically (`autostash`). See [this blog post](https://cscheng.info/2017/01/26/git-tip-autostash-with-git-pull-rebase.html) on more information about the `git rebase --autostash`-functionality.

```bash
git config rebase.autoStash true
```

Create a second remote called `upstream` that points at the official OGS repository and fetch from it:

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

You will also need to install `clang-format`:

<div class='win'>

Install clang (which contains `clang-format`) with the [official installer](https://prereleases.llvm.org/win-snapshots/LLVM-12.0.0-6923b0a7-win64.exe)

</div>

<div class='linux'>

```bash
sudo apt-install clang-format
```

</div>

<div class='mac'>

```bash
brew install clang-format
```

</div>

## Step: Working on a new feature

You only have to follow the above steps once. From then on, whenever you want to work on a new feature, you can more easily interact with the remote repositories.

Make sure that your local repository is up-to-date with the `upstream` repository:

```bash
git fetch upstream
```

Create a branch `feature-name` off of `upstream/master`-branch to work on a new feature, and check out the branch:

```bash
git checkout -b feature-name upstream/master
```

```mermaid
gitGraph
  commit
  commit
  branch feature-name
  commit
```

----

To keep up to date with the developments in the official repository it is recommended to rebase your feature-branch regularly (at least weekly) with:

```bash
git fetch upstream
git rebase upstream/master
```

<div class="note">

This can potentially lead to conflicts, which [have to be resolved](https://www.git-tower.com/learn/git/faq/solve-merge-conflicts).

</div>

----

Now after you implemented the feature and committed your work you can push the new commits to the `feature-name`-branch on your GitLab fork:

```bash
git push -u origin feature-name  # -u is required only first time to set up the remote-tracking.
```

<div class="note">

In case you already have pushed your branch before and also rebased on `master` afterwards you may [need to do a force `push`](https://www.git-tower.com/learn/git/faq/git-force-push)!

```bash
git push --force origin feature-name
```

</div>

If your work is done submit a [merge request](https://gitlab.opengeosys.org/ogs/ogs/-/merge_requests/new).

----

Again this triangular workflow is summarized with this picture:
![Workflow visualization](https://cloud.githubusercontent.com/assets/1319791/8943755/5dcdcae4-354a-11e5-9f82-915914fad4f7.png)

<!-- vale off -->

> […] the Forking Workflow requires two remotes—one for the official repository, and one for the developer’s personal server-side repository. While you can call these remotes anything you want, a common convention is to use origin as the remote for your forked repository […] and upstream for the official repository.
>
> <cite><a href="https://www.atlassian.com/git/tutorials/comparing-workflows/forking-workflow">www.atlassian.com/git/tutorials/comparing-workflows/forking-workflow</a> </cite>

<!-- vale on -->
