+++
date = "2018-02-26T11:00:13+01:00"
title = "Set Up GitLab"
author = "Lars Bilke"
weight = 1002

[menu]
  [menu.devguide]
    parent = "getting-started"
+++

## Introduction

[GitLab](https://gitlab.com) is a web-based development and collaboration tool similar to [GitHub](https://github.com). We self-host GitLab at <https://gitlab.opengeosys.org> and migrated our development workflows from GitHub in June 2020. Development takes place in the [ogs-group](https://gitlab.opengeosys.org/ogs) and the authoritative repository is at [ogs/ogs](https://gitlab.opengeosys.org/ogs/ogs).

## Setup an account

- Creating a GitLab account can be done by simply using your existing GitHub account: click the GitHub logo (octocat) on the [Gitlab sign-in page](https://gitlab.opengeosys.org/users/sign_in)
  ![GitLab login page](../gitlab-login.png)
  - You will be redirected to GitHub (please login there) and asked for authorization.
  - Your new user account will be blocked at first, please let us know we will unblock it

## Setup a password for cloning over https

To clone a repository via `https://`-protocol you need to [set up an account password](https://gitlab.opengeosys.org/profile/password/edit) and use this during `git clone https://...`.

## Upload a SSH key for cloning over SSH

To clone a repository via SSH-protocol (`git clone git@gitlab...`) you need to [upload your SSH public key](https://gitlab.opengeosys.org/profile/keys).
