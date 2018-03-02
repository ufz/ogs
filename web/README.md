## The big picture

- Development related content such as developer guide, benchmark documentation, tools description, ... will be simple Markdown files in e.g. [/web/content/docs](https://github.com/bilke/ogs/blob/web-hugo/web/content/docs/benchmarks/elliptic/groundwater-flow-neumann.md)
- You can preview documentation locally with [Hugo](https://gohugo.io) – a static site generator
- Other content such as news, blog posts, articles can be authored in a graphical environment and is not part of the code but hosted with a content management system (currently a prototypical implementation uses [Contenful](https://www.contentful.com/) but maybe we can even stick to the current system [Craft CMS](https://craftcms.com/) once version 3 is released). This will lower barriers for non-technical people to contribute to the site. The content gets imported automatically on Jenkins or locally (see [web/README.md](https://github.com/bilke/ogs/blob/web-hugo/web/README.md) for details)
- You can [mark](https://github.com/bilke/ogs/blob/web-hugo/ProcessLib/GroundwaterFlow/CMakeLists.txt#L80) benchmarks to be automatically and interactively visualized [in the documentation](https://github.com/bilke/ogs/commit/d4fc7d94a3821a6b4483a1d7aeaabd6ee391c449#diff-2f5b1ac2a759aa09b2d3f5cc1ece45ceR108) inside your [browser](https://dev.opengeosys.org/docs/benchmarks/elliptic/groundwater-flow-neumann/#results-and-evaluation) via [vtk.js](https://kitware.github.io/vtk-js/)! 🍻 CURRENTLY DISABLED!

## Requirements

- Download [Hugo](https://gohugo.io/#action) and put it in your `PATH`
- [Install Pandoc](https://pandoc.org/installing.html)

## Getting started

- Inside the source-directory `ogs/web`:
  - Run `hugo server`
  - Open http://localhost:1313

As you make modifications to the site it will be rebuild and the page in the browser gets reloaded.

## How-Tos

### Create a new page

By using `hugo new` you can create a new page with the correct frontmatter for that kind of page:

```bash
hugo new docs/benchmarks/elliptic/groundwater-flow-dirichlet.pandoc
```

- `--kind` is one the `archetypes/*`
- path is relative to `content/` and determines the URL of the page

Or you can simply create a new `.pandoc`-file in the correct location and fill it by yourself.

### Setup navigation for a page

The are submenus (shown in the left sidebar) for specific sections such as for benchmarks. The submenus consist of groups (e.g. *Elliptic*) and page entries. Groups are defined in `config.toml`:

```toml
[[menu.benchmarks]]
  name = "Elliptic"
  identifier = "elliptic"
  weight = 1
```

To add your page to a group as an entry add the following frontmatter:

```toml
weight = 101

[menu]
  [menu.benchmarks]
    parent = "elliptic"
```

`weight` specifies the order of groups and pages in ascending order (top -> down).

### Write a page

We use [Pandoc Markdown](https://pandoc.org/MANUAL.html#pandocs-markdown) for the actual content.

It is an enhanced version of the [original Markdown](http://daringfireball.net/projects/markdown/). Please consult both guides!

#### Images

Use regular Markdown syntax:

```md
![](../square_1e2_neumann_gradients.png)
```

The path to the image can be absolute (by preceding with `/`) or relative. The relative path starts at your current URL. If your image is in the same directory as your `.pandoc`-file you have to prefix your path with `../` as in the example above.

See the [Pandoc Help](https://pandoc.org/MANUAL.html#images) for more options on e.g. image size and captions.

#### Visualizations

CURRENTLY DISABLED!

Use shortcode `vis`:

```go
{{< vis path="Elliptic/square_1x1_GroundWaterFlow/square_1e2_pcs_0_ts_1_t_1.000000.vtu" [height="300"] >}}
```

`path` is relative to `ogs/web/static/vis/` and points to a converted data set, see below.

Optional parameters:

- `height` - In px

----

You can convert VTK output files to vtk.js format with:

```bash
npm run convert -- -e -i vtk-file.vtu -o output/path
```

Run `npm run convert` to get a list of possible arguments.

#### Equations

Equations can be set with typical LaTeX syntax by using [MathJax](https://www.mathjax.org/). Blocks are defined by `$$` at the beginning and `$$` at the end of the block. Inline math uses one `$` as the delimiter. For more usage instructions see the [MathJax LaTeX help](http://docs.mathjax.org/en/latest/tex.html).

#### Bibliography references

Bibliography items from *Documentation/bibliography.bib* can be referenced by their id with the `bib`-shortcode:

```go
{{< bib id="Kolditz2012" >}}
```

The bib-file has to be converted into a json-file with the [pandoc-citeproc](https://github.com/jgm/pandoc-citeproc)-tool:

```sh
[cd to ogs/web]
pandoc-citeproc --bib2json ../Documentation/bibliography.bib > data/bibliography.json
```

This json-file is then used by the shortcode.

---

## Advanced topics

### Optional requirements

- Install [Yarn](https://yarnpkg.com/en/docs/install); for downloading required JavaScript & CSS development packages
- Install Python and `pip`; for getting content from [Contentful](https://contentful.com)
- Install [ParaView](http://www.paraview.org/download/); for converting VTK output files to vtk.js-format for interactive web visualization, check if `pvpython` is either in the `PATH` or ParaView is installed in `/Applications/` on macOS or `/usr/local/opt/paraview` on Linux. CURRENTLY DISABLED!

### CSS & JavaScript development

- Inside the source-directory `ogs/web`:
  - Install packages with `yarn`
  - Run `webpack --watch`

    - Re-run CMake and build the `ctest`-target, *OPTIONAL* for benchmark visualizations
  - Run `npm run build` to build the site which is created in `public/`

### Importing CMS content (Optional)

- Install Python packages with `pip install -r requirements.txt`
- In `ogs/web/import` rename `secret_example.py` to `secret.py` and fill in `accessToken`.
- In `ogs/web/import` run `python import.py`

This fetches articles from the CMS to e.g. `ogs/web/data/news.json`.

### Update search index

```bash
ALGOLIA_WRITE_KEY=XXX node_modules/.bin/hugo-algolia --toml -s
```

### Used components

- [Hugo](https://gothugo.com) - Web site generator
- [Contenful](https://www.contentful.com/) -  API-based CMS for news, articles, ..
- [Tailwind](https://tailwindcss.com/docs/what-is-tailwind) - CSS framework
- [vtk.js](https://kitware.github.io/vtk-js/) - 3D Visualizations, CURRENTLY DISABLED!
- [webpack](https://webpack.github.io/) - Packaging JavaScript & CSS
- [FontAwesome](https://fontawesome.com) - Icons, see [icon search](https://fontawesome.com/icons?d=gallery)
