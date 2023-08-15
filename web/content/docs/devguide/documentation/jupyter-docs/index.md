+++
date = "2022-06-28T15:28:13+01:00"
title = "Jupyter notebooks as documentation and tests"
author = "Lars Bilke"
weight = 1025

[menu.devguide]
parent = "development-workflows"
+++

## The big picture

[Jupyter notebooks](https://jupyter.org) are interactive computing environments where prose and code can be combined. In the OGS project notebooks can be used to define complex benchmark workflows and its results can be converted to be shown on the OGS web page (see an [example here](/docs/benchmarks/small-deformations/linear_disc_with_hole/)).

## Create a new notebook

Create a new notebook file in `Tests/Data` (either as regular `.ipynb`-files or as Markdown-based notebooks via [Jupytext](https://jupytext.readthedocs.io/en/latest)). See examples:

- [SimpleMechanics.ipynb](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Mechanics/Linear/SimpleMechanics.ipynb) (regular `.ipynb`-notebook)
- [Linear_Disc_with_hole.md](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Mechanics/Linear/DiscWithHole/Linear_Disc_with_hole.md) (Jupytext-based notebook in Markdown)

## Add web meta information

If the notebook result should appear as a page on the web documentation a frontmatter with some meta information (similar to [regular web pages]({{< ref "web-docs.md" >}})) is required as the first cell in the notebook:

- Add a new cell and move it to the first position in the notebook
- Change the cell type to `raw`!
- Add meta information, conclude with a end-of-file marker (`<!--eofm-->`) e.g.:

  ```md
  title = "SimplePETSc"
  date = "2021-11-09"
  author = "Lars Bilke"
  image = "optional_gallery_image.png"
  web_subsection = "small-deformations"
  <!--eofm-->
  ```

<div class='note'>

In Jupytext-based notebooks you can add the frontmatter within the `<!-- #raw -->`- and `<!-- #endraw -->`-markers:

```md
<!-- #raw -->
title = "Linear elasticity: disc with hole"
date = "2022-04-27"
author = "Linda Günther, Sophia Einspänner, Robert Habel, Christoph Lehmann and Thomas Nagel"
web_subsection = "small-deformations"
<!--eofm-->
<!-- #endraw -->
```

`web_subsection` needs to be set to a sub-folder in [web/content/docs/benchmarks](https://gitlab.opengeosys.org/ogs/ogs/-/tree/master/web/content/docs/benchmarks) (if not set the notebook page will not be linked from navigation bar / benchmark gallery on the web page).

</div>

## Notebook setup

### Python cells

- Do not use machine-specific or absolute paths! See the following example to set up notebook output paths:

  ```python
  import os

  # On CI out_dir is set to the notebooks directory inside the build directory
  # similar to regular benchmark tests. On local testing it will output to the
  # notebooks source directory under a _out-subdirectory.
  out_dir = os.environ.get('OGS_TESTRUNNER_OUT_DIR', '_out')
  if not os.path.exists(out_dir):
      os.makedirs(out_dir)

  # ...
  # Run ogs; get input data from current directory; write to `out_dir`
  ! ogs my_project.prj -o {out_dir} > {out_dir}/log.txt

  # OR with ogs6py:
  # ... setup model ...
  model.run_model(logfile=os.path.join(out_dir, "log.txt"), args=f"-o {out_dir}")

  # Verify results; on failure assert with:
  assert False
  # or
  raise SystemExit()
  ```

- Do not write anything into the source directories. Use an `out_dir` as above.
- Assume that `ogs` and other tools are in the `PATH`.

### Markdown cells

- Do not use HTML inside Markdown blocks.
- Static images e.g. for the gallery image or to be used in Markdown cells have to be located in either `images`- or `figures`-subdirectories beneath the Notebook file! Otherwise they will not show up in the web site.
  - For image captions add a title in quotation marks after the image path, e.g.

    ```md
    ![Alt text](figures/my_image.png "This will be the image caption.")
    ```

  - Please note that in merge request web previews static images are not shown at all.

## Execution environment

In CI the notebooks are executed with all dependencies installed into a virtual environment in the build directory. The installed packages are defined in `Test/Data/requirements.txt`. The same setup can be enabled locally by setting the CMake option `OGS_USE_PIP=ON`. E.g.

```bash
cmake --preset release -DOGS_USE_PIP=ON    # Creates the virtual environment
source ../build/release/.venv/bin/activate # Activates the virtual environment
jupyter lab                                # Starts Jupyter Lab
```

## Register with CTest

Add the notebook to CTest with e.g.:

```cmake
if(NOT OGS_USE_PETSC)
    NotebookTest(NOTEBOOKFILE Mechanics/Linear/SimpleMechanics.ipynb RUNTIME 10)
endif()
```

Then e.g. run all notebook test (`-R nb`) in parallel (`-j 4`) with:

```bash
# cd into build directory
source .venv/bin/activate # Is created with OGS_USE_PIP=ON, see above note on environment.
ctest -R nb -j 4 --output-on-failure
```

## Advanced topics

### Jupytext usage

If you use the [execution environment](#execution-environment) [Jupytext](https://jupytext.readthedocs.io/en/latest) is already configured and its usage is transparent:

- Double-click on a markdown file will open it as a Notebook
- Upon saving or executing a linked `.ipynb`-file is created in the background which stores outputs
- You still edit the Markdown file but don't notice the difference to regular notebooks in the Lab UI

### Run a notebook in BinderHub

On the web site or MR web previews on pages generated by a notebook there is a new banner:

![Notebook web banner with BinderHub launch button](binderhub-button.png)

- Click the button to launch the notebook in BinderHub.
- The environment running in BinderHub is defined in [`bilke/binder-ogs-requirements` at GitHub](https://github.com/bilke/binder-ogs-requirements)
- When clicking the link it launches a Jupyter Lab instance pre-configured with ogs [via wheel](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/requirements-ogs.txt#L2), clones the current ogs repo in it and opens the respective notebook ready to run. Please note that startup times may be several minutes and the computing resources are limited (1 core, 2GB RAM). For improved performance we would need to setup own infrastructure. Also currently only works for serial ogs configurations.

### PyVista notebooks on headless Linux systems

PyVista (or VTK) requires a windowing environment for rendering. You can provide a virtual window with `xvfb-run`:

```bash
sudo apt install libgl1-mesa-glx xvbf # install xvfb

xvfb-run -a ctest [...] # provide a virtual window to the ctest-run
```
