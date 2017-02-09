## Requirements

- Download [Hugo](https://gohugo.io/#action) and put it in your `PATH`
- Install [Node.js](https://nodejs.org/en/), check if `npm` is in your `PATH`
- Install [Yarn](https://yarnpkg.com/en/docs/install), *OPTIONAL* faster alternative to `npm install`
- Install Python and `pip`, *OPTIONAL* for getting content from [Contentful](https://contentful.com)
- Install [ParaView](http://www.paraview.org/download/), *OPTIONAL* for converting VTK output files to vtk.js-format for interactive web visualization, check if `pvpython` is either in the `PATH` or ParaView is installed in `/Applications/` on macOS or `/usr/local/opt/paraview` on Linux

## Getting started

- Inside the source-directory `ogs/web`:
  - Install Node packages with `npm install` (or via `yarn`)
  - Install Python packages with `pip install -r requirements.txt`, *OPTIONAL* for getting content from [Contentful](https://contentful.com)
  - Run `npm run build` to build the site which is created in `public/`

Alternatively you can simply run `make web` inside your build directory to install everything and build the site (this requires that `python`, `pip` and `npm` or `yarn` was found by CMake).

## Develop site with live-preview

In `ogs/web` run

    hugo server

In your browser go to http://localhost:1313. As you make modifications to the site it will be rebuild and the page in the browser gets reloaded.

If you want to modify css or javascript run `gulp` in another terminal:

    npm run gulp

## Importing CMS content (Optional)

In `ogs/web/import` rename `secret_example.py` to `secret.py` and fill in `accessToken`.

In `ogs/web/import` run

    python import.py

This fetches articles from the CMS to e.g. `ogs/web/data/news.json`.

## Used components

- [Hugo](https://gothugo.com) - Static site generator for technical documentation
- [Contenful](https://www.contentful.com/) -  API-based CMS for news, articles, ..
- [flexboxgrid](http://flexboxgrid.com/) - CSS grid
- [vtk.js](https://kitware.github.io/vtk-js/) - 3D Visualizations
- [webpack](https://webpack.github.io/) - Packaging JavaScript
- [gulp](http://gulpjs.com/) - Automation toolkit

## How-Tos

### Images

Use shortcode `img`:

```go
{{< img src="../square_1e2_neumann_gradients.png" >}}
```

`src` can be absolute (by preceding with `/`) or relative. The relative path starts at your current URL. If your image is in the same directory as your `.md`-file you have to prefix your path with `../` as in the example above.

Optional parameters:

- `title` - Image caption
- `class` - CSS class
- `alt` - Alt text

### Visualizations

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
