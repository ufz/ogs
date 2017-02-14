+++
title = "About"
author = "Lars Bilke"
date = "2017-01-13T14:24:23+01:00"

[menu]
    [menu.main]
        name = "About"
        url = "/about/"
        weight = 2
+++

See [web/README.md](https://github.com/ufz/ogs/blob/master/web/README.md) for getting started on developing this site!

## The big picture

- Development related content such as developer guide, benchmark documentation, tools description, ... will be simple Markdown files in e.g. [/web/content/docs](https://github.com/bilke/ogs/blob/web-hugo/web/content/docs/benchmarks/elliptic/groundwater-flow-neumann.md)
- You can preview documentation locally with [Hugo](https://gohugo.io) – a static site generator
- Other content such as news, blog posts, articles can be authored in a graphical environment and is not part of the code but hosted with a content management system (currently a prototypical implementation uses [Contenful](https://www.contentful.com/) but maybe we can even stick to the current system [Craft CMS](https://craftcms.com/) once version 3 is released). This will lower barriers for non-technical people to contribute to the site. The content gets imported automatically on Jenkins or locally (see [web/README.md](https://github.com/bilke/ogs/blob/web-hugo/web/README.md) for details)
- You can [mark](https://github.com/bilke/ogs/blob/web-hugo/ProcessLib/GroundwaterFlow/CMakeLists.txt#L80) benchmarks to be automatically and interactively visualized [in the documentation](https://github.com/bilke/ogs/commit/d4fc7d94a3821a6b4483a1d7aeaabd6ee391c449#diff-2f5b1ac2a759aa09b2d3f5cc1ece45ceR108) inside your [browser](https://dev.opengeosys.org/docs/benchmarks/elliptic/groundwater-flow-neumann/#results-and-evaluation) via [vtk.js](https://kitware.github.io/vtk-js/)! 🍻

### Notable changes

- [git-lfs](https://git-lfs.github.com/) is used for image files in the documentation section, so this is a new requirement if you want to edit the docs
- New build target `web` can be used to build the docs
- Site is built by Jenkins on every build (also on PRs), just click in the link called [Web](https://jenkins.opengeosys.org/job/User/job/bilke/job/ogs/job/web-hugo/Web/) on the sidebar
- Site built from `ufz/master` is served at https://dev.opengeosys.org at the moment

## Requirements

All of the following should be optional if you don't want to deal with the site.

- [Hugo](https://gohugo.io), [NodeJS](https://nodejs.org/en/) for web site generation
- Python, pip and super secret credentials for importing content from CMS
- git-lfs to checkout images
- ParaView installed for converting benchmark result files to vtk-js for web visualization
