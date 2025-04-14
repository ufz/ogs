+++
date = "2021-09-09T15:00:13+01:00"
title = "Jupyter Notebooks"
author = "Lars Bilke"
weight = 61
+++

<!-- TODO: Consider to move this section out of **BASICS** and to devote an extra first-order section to the Python-bindings of OGS and how to operate OGS via Python and Jupyter Notebooks.-->

[Jupyter Notebooks](https://jupyter.org) are documents which can contain live (Python) code, equations, visualizations and
narrative text and can be used as an intuitive interface for OGS projects. The following video gives an introduction to using
OpenGeoSys with Jupyter Notebooks:

{{< youtube eihNKjK-I-s >}}

## Jupyter Notebooks container environments

You can start with pre-defined container environment from one of the images

### Usage

With [Docker]({{< ref "container.md#with-docker" >}}):

```bash
docker run --rm -p 8888:8888 -v $PWD:/home/jovyan/work --user `id -u $USER` \
    --group-add users quay.io/jupyter/scipy-notebook
```

This mounts your current directory into `~/work` inside the container.

<div class="note">

#### <i class="fab fa-windows"></i> Windows notes

The above command only works when you run Docker from within a WSL2 Linux shell!

- Install [Docker Desktop on Windows](https://docs.docker.com/desktop/windows/install/), it automatically configures WSL2.
- Install a Linux distribution from the Microsoft App Store. We recommend [Ubuntu 22.04](https://apps.microsoft.com/search?query=ubuntu).
- In the Docker Desktop application under *Settings / Resources / WSL integration* add your Linux distribution.
- Open a command prompt in your Linux distribution (At the start menu type the name of the distribution) and run the container.
  - If your current working contains spaces write out `$PWD`, e.g.:

    ```bash
    ... -v /c/Users/My\ Name/working/directory:/home/jovyan/work ...
    ```

</div>

---

With [Singularity]({{< ref "container.md#with-singularity" >}}):

```bash
singularity run docker://quay.io/jupyter/scipy-notebook
```

Open the specified URL shown in the command output in your browser, e.g.

```bash
http://127.0.0.1:8888/lab?token=xxx
```

You may have to modify the IP address if this is running on a remote machine.

### Install ogs

Click on *Terminal* and install OGSTools (which also installs ogs itself):

```bash
pip install ogstools[ogs]
```

You install other Python packages too. Please note that this is a temporary installation. If you stop the container the manually installed packages will be deleted.
