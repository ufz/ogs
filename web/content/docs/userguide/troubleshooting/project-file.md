+++
date = "2018-11-14T11:00:13+01`:00"
title = "Project file syntax"
author = "Lars Bilke"
weight = 3

[menu]
  [menu.userguide]
    parent = "troubleshooting"
+++

## Check project file syntax with `xmllint`

Project files `.prj` have to be valid XML documents. You can check the formatting with the [`xmllint`-tool](http://xmlsoft.org/xmllint.html):

```bash
xmllint --noout myproj.prj
```

### Install `xmllint`

<div class='win'>

We recommend to install via [Chocolatey](https://chocolatey.org):

```powershell
choco install xsltproc
```

<div class='note'>

### <i class="far fa-info-circle"></i> Alternative installation

Another method is to use the [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10) where you can simply install Linux packages:

```bash
sudo apt-get install libxml2-utils
```

</div>

</div>

<div class='linux'>

```bash
sudo apt-get install libxml2-utils
```

</div>

<div class='mac'>
`xmllint` is part of the Homebrew xmlstarlet package:

```bash
brew install xmlstarlet
```

</div>
