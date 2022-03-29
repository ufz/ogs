+++
date = "2018-02-23T15:28:13+01:00"
title = "Get the source code"
author = "Lars Bilke"
weight = 1003

[menu]
  [menu.devguide]
    parent = "getting-started"
+++

<div class='note'>

<i class="far fa-exclamation-triangle"></i> This page describes how to get the source code as a simple download or git clone. The information on the forking workflow with `git` has been moved to the [Development workflows]({{< ref "setup-fork.md" >}})-section.

</div>

## Option: Download the source code

To get the latest source code simply download it from [repository website](https://gitlab.opengeosys.org/ogs/ogs) and extract the archive:

![](zip-download.png)

You can also download and extract with the command line:

```bash
wget https://gitlab.opengeosys.org/ogs/ogs/-/archive/master/ogs-master.tar.gz
tar xf ogs-master.tar.gz
```
----
## Option: Clone the source code repository with Git

First you need to get the clone url:

![](git-url.png)

Then clone the repository with `git`:

```bash
git clone https://gitlab.opengeosys.org/ogs/ogs.git
```
