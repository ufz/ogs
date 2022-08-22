+++
date = "2022-06-28T15:28:13+01:00"
title = "Jupyter notebooks as documentation"
author = "Lars Bilke"
weight = 1025

[menu.devguide]
parent = "development-workflows"
+++

## The big picture

Jupyter notebooks in `Tests/Data` are automatically executed and converted to web pages in the benchmark documentation section.

## Add web meta information

Similar to regular web documentation pages the notebook requires to have a frontmatter with some meta information as the first cell in the notebook:

- Add a new cell and move it to the first position in the notebook
- Change the cell type to `raw`!
- Add meta information, conclude with a end-of-file marker (`<!--eofm-->`) e.g.:

  ```toml
  title = "SimplePETSc"
  date = "2021-11-09"
  author = "Lars Bilke"
  notebook = "Notebooks/SimplePETSc.ipynb"
  image = "optional_gallery_image.png"
  web_subsection = "elliptic"
  <!--eofm-->
  ```

Static images e.g. for the gallery image or to be used in Markdown cells have to located in either `images`- or `figures`-subdirectories beneath the Notebook file! Otherwise they will not show up in the web site.

## General advice

### Python cells

- See [Notebook testing]({{< relref "test-data#notebook-testing" >}}) for details on how to setup ogs execution (especially input data path handling) inside the notebook. See also the [SimplePETSc.ipynb](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/Notebooks/SimplePETSc.ipynb)-notebook as an example.
- Do not use machine-specific or absolute paths!
- Assume that ogs and other tools are in the `PATH`.

### Markdown cells

- Do not use HTML inside Markdown blocks.
- For image captions add a title in quotation marks after the image path, e.g.

  ```md
  ![Alt text](figures/my_image.png "This will be the image caption.")
  ```

- Please use `images`- or `figures`-subdirectories!
- Please note that in merge request web previews images are not shown at all.
