## Getting started

- Download [Hugo](https://gohugo.io/#action) and put it in your `PATH`
- Install gulp-cli globally with `npm install --global gulp-cli`, *OPTIONAL* for SCSS and Javascript
- Install Node packages with `npm install`, *OPTIONAL* for SCSS and Javascript
- Install Python packages with `pip install -r requirements.txt`, *OPTIONAL* for getting content from [Contentful](https://app.contentful.com/spaces/4nuqzxntzxks)

## Start servers and watchers

    hugo server

If you want to modify css or javascript run `gulp` in another terminal:

    gulp

## Build site

    (cd import; python import.py) # Optional for fetching content from Contentful
    gulp build
    hugo

Test by locally serving via [Caddy](https://caddyserver.com):

    caddy

## Used components

- [Hugo](https://gothugo.com) - Static site generator for technical documentation
- [Contenful](https://www.contentful.com/) -  API-based CMS for news, articles, ..
- [flexboxgrid](http://flexboxgrid.com/) - CSS grid
- [vtk-js](https://kitware.github.io/vtk-js/) - 3D Visualizations
- [webpack](https://webpack.github.io/) - Packaging JavaScript
- [gulp](http://gulpjs.com/) - Automation toolkit

## Dump

### Serve converted meshes with S3

- Set CORS: http://stackoverflow.com/a/19557452/80480
- Upload to e.g. http://s3-eu-central-1.amazonaws.com/jenkins-dump
- Eventually add `crossdomain.xml` to the root of the bucket:

```xml
<?xml version="1.0"?>
<!DOCTYPE cross-domain-policy SYSTEM
"http://www.macromedia.com/xml/dtds/cross-domain-policy.dtd">
<cross-domain-policy>
  <allow-access-from domain="*" secure="false" />
</cross-domain-policy>
```
