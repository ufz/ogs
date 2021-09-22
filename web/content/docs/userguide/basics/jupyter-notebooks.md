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

You can use a pre-defined container environment which currently contains:

- The Jupyter Notebook / Lab application
- The latest OpenGeoSys application (with MFront support) and tools
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
- Jupyter-related tools:
  - [nbconvert](https://nbconvert.readthedocs.io) — Format conversion
  - [nbdime](https://nbdime.readthedocs.io) — Diffs for notebooks
  - [nb2hugo](https://github.com/bilke/nb2hugo/tree/ogs) — Notebook to www.opengeosys.org markdown
### Usage

With [Docker]({{< ref "container.md#with-docker" >}}):

```bash
docker run --rm -p 8888:8888 -v $PWD:/home/jovyan/work --user `id -u $USER` \
    --group-add users registry.opengeosys.org/ogs/ogs/ogs-serial-jupyter
```

This mounts your current directory into `~/work` inside the container.

---

With [Singularity]({{< ref "container.md#with-singularity" >}}):

```bash
singularity run docker://registry.opengeosys.org/ogs/ogs/ogs-serial-jupyter
```

Open the specified URL shown in the command output in your browser, e.g.

```
http://127.0.0.1:8888/lab?token=xxx
```

You may have to modify the IP address if this is running on a remote machine.

### Adding additional Python packages

In a running container you can install additional Python packages with the Jupyter [magic command `%pip`](https://ipython.readthedocs.io/en/stable/interactive/magics.html#magic-pip):

```python
%pip install [package name]
```

Please note that this is a temporary installation. If you stop the container the environment is destroyed.

### Rendering with PyVista

When using [PyVista][pyvista] the container uses the (interactive) [pythreejs](https://docs.pyvista.org/user-guide/jupyter/pythreejs.html) rendering backend per default. If you want to output static images there are a couple of ways to configure the 'static' rendering backend:

- Globally as an environment variable: `PYVISTA_JUPYTER_BACKEND=static`. This can be defined when starting the container:
  ```bash
  docker run ... -e PYVISTA_JUPYTER_BACKEND=static ...
  ```
- Globally as a pyvista setting:
  ```python
  import pyvista as pv
  pv.set_jupyter_backend('static')
  ```
- Locally on each plot:
  ```python
  import pyvista as pv
  ...
  pv.show(jupyter_backend='static')
  ```

[pyvista]: https://docs.pyvista.org
