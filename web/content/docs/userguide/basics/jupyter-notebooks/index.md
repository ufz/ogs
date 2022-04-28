+++
date = "2021-09-09T15:00:13+01:00"
title = "Jupyter Notebooks"
author = "Lars Bilke"
weight = 5

[menu]
  [menu.userguide]
    parent = "basics"
+++

[Jupyter Notebooks](https://jupyter.org) are documents which can contain live (Python) code, equations, visualizations and narrative text and can be used as an intuitive interface for OGS projects. The following video gives an introduction to using OpenGeoSys with Jupyter Notebooks:

{{< youtube eihNKjK-I-s >}}

## Jupyter Notebooks container environments

You can use a pre-defined container environment.

Image `registry.opengeosys.org/ogs/ogs/ogs-serial-jupyter` contains:

- The Jupyter Notebook / Lab application
- The latest OpenGeoSys application with MFront-support and tools
- A set of Python packages:
  - [ogs6py](https://github.com/joergbuchwald/ogs6py) — OGS model manipulation
  - [VTUInterface](https://github.com/joergbuchwald/VTUinterface) — VTU / PVD IO
  - [h5py](https://docs.h5py.org/en/latest/index.html) — HDF5 IO
  - [MFront python bindings](http://tfel.sourceforge.net/mfront-python.html) – Material model manipulation
  - [matplotlib](https://matplotlib.org) — Plotting
  - [numpy](https://numpy.org) — Scientific computing
  - [pandas](https://pandas.pydata.org) — Data analysis
  - [scipy](https://docs.scipy.org/doc/scipy/reference/) — Scientific computing
  - [vtk](https://pypi.org/project/vtk/) — Visualization
  - [PyVista][pyvista] — Visualization
  - [Snakemake](https://snakemake.github.io) — Workflow management
- Jupyter-related tools:
  - [nbconvert](https://nbconvert.readthedocs.io) — Format conversion
  - [nbdime](https://nbdime.readthedocs.io) — Diffs for notebooks
  - [nb2hugo](https://github.com/bilke/nb2hugo/tree/ogs) — Notebook to www.opengeosys.org markdown

Image `registry.opengeosys.org/ogs/ogs/ogs-petsc-jupyter` additionally contains:

- PETSc-support

### Usage

With [Docker]({{< ref "container.md#with-docker" >}}):

```bash
docker run --rm -p 8888:8888 -v $PWD:/home/jovyan/work --user `id -u $USER` \
    --group-add users registry.opengeosys.org/ogs/ogs/ogs-serial-jupyter
```

This mounts your current directory into `~/work` inside the container. Use image `registry.opengeosys.org/ogs/ogs/ogs-petsc-jupyter` for PETSc-support!


<div class="note">

#### <i class="fab fa-windows"></i> Windows notes

The above command only works when you run Docker from within a WSL2 Linux shell!

- Install [Docker Desktop on Windows](https://docs.docker.com/desktop/windows/install/), it automatically configures WSL2.
- Install a Linux distribution from the Microsoft App Store. We recommend [Ubuntu 20.04](https://www.microsoft.com/en-us/p/ubuntu-2004-lts/9n6svws3rx71).
- In the Docker Desktop application under *Settings / Resources / WSL integration* add your Linux distribution.
- Open a command prompt in your Linux distribution (At the start menu type the name of the distribution) and run the container.
    - If your current working contains spaces write out `$PWD`, e.g.:
    ```
    ... -v /c/Users/My\ Name/working/directory:/home/jovyan/work ...
    ```

</div>

---

With [Singularity]({{< ref "container.md#with-singularity" >}}):

```bash
singularity run docker://registry.opengeosys.org/ogs/ogs/ogs-serial-jupyter
```

Open the specified URL shown in the command output in your browser, e.g.

```bash
http://127.0.0.1:8888/lab?token=xxx
```

You may have to modify the IP address if this is running on a remote machine.

<div class="note">

#### <i class="fab fa-windows"></i> Specific OGS version

You can append a version number to the image name (applies both to Docker and Singularity) to get an image for a specific OGS release (starting with 6.4.1):

```
singularity run docker://registry.opengeosys.org/ogs/ogs/ogs-serial-jupyter:6.4.1
```

Available images are [listed on GitLab](https://gitlab.opengeosys.org/ogs/ogs/container_registry/79).

</div>

### Browsing notebooks on GitLab

In the file browser on the left-hand side of the Jupyter Lab interface there is a GitLab-tab which allows for browsing and opening notebooks from the [ogs/ogs](https://gitlab.opengeosys.org/ogs/ogs)-repository. You can directly modify and execute a notebook, but the notebook is not saved back to GitLab. You can change the browsed repository by typing into the top text field.

If you would like to us this with private repositories you have to supply an [access token](https://gitlab.opengeosys.org/-/profile/personal_access_tokens) at container start-up:

```bash
docker run --rm -p 8888:8888 -v $PWD:/home/jovyan/work --user `id -u $USER` \
    --group-add users registry.opengeosys.org/ogs/ogs/ogs-serial-jupyter \
    --GitLabConfig.access_token="< YOUR_ACCESS_TOKEN >"
```

### Adding additional Python packages

In a running container you can install additional Python packages with the Jupyter [magic command `%pip`](https://ipython.readthedocs.io/en/stable/interactive/magics.html#magic-pip):

```python
%pip install [package name]
```

Please note that this is a temporary installation. If you stop the container the environment is destroyed.

### Rendering with PyVista

When using [PyVista][pyvista] the container uses the (interactive local rendering) [pythreejs](https://docs.pyvista.org/user-guide/jupyter/pythreejs.html) rendering backend per default. Other backends are currently not supported.

[pyvista]: https://docs.pyvista.org
