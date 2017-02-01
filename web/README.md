## Start servers

    hugo server

If you want to modify css or javascript run `gulp` in another terminal:

    gulp
    
## Build site

    (cd import; python import.py)
    gulp build
    hugo

Test by locally serving via [Caddy](https://caddyserver.com):

    caddy

## Used components

- [Hugo](https://gothugo.com)
- [flexboxgrid](http://flexboxgrid.com/)


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
