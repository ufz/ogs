+++
date = "2018-11-14T15:00:13+01:00"
title = "Jupyter Notebooks in a container"
author = "Lars Bilke"
weight = 5

[menu]
  [menu.userguide]
    parent = "basics"
+++

[Jupyter Notebooks](https://jupyter.org) are documents which can contain live (Python) code, equations, visualizations and narrative text and can be used as an intuitive interface for OGS projects.

{{< youtube eihNKjK-I-s >}}

You can use a pre-defined container environment which currently contains:

- The Jupyter Notebook application
- The latest OpenGeoSys application and tools
- A set of Python packages:
  - [ogs6py](https://github.com/joergbuchwald/ogs6py)
  - [VTUInterface](https://github.com/joergbuchwald/VTUinterface)
  - [matplotlib](https://matplotlib.org)
  - [numpy](https://numpy.org)
  - [pandas](https://pandas.pydata.org)
  - [scipy](https://docs.scipy.org/doc/scipy/reference/)
  - [vtk](https://pypi.org/project/vtk/)

## Usage

With [Docker]({{< ref "container.md#with-docker" >}}):

```bash
docker run --rm -p 8888:8888 -v $PWD:/home/jovyan/work --user `id -u $USER` \
    --group-add users registry.opengeosys.org/ogs/ogs/ogs-serial-jupyter
```

With [Singularity]({{< ref "container.md#with-singularity" >}}):

```bash
singularity run docker://registry.opengeosys.org/ogs/ogs/ogs-serial-jupyter
```

Open the specified URL shown in the command output in your browser.
