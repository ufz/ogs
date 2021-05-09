+++
date = "2018-11-14T11:00:13+01`:00"
title = "General"
author = "Lars Bilke"
weight = 22
toc = true

[menu]
  [menu.userguide]
    parent = "troubleshooting"
+++

## Data Explorer

### XSDError: Loaded schema file is invalid

You may encountering the following error (or similar) on opening `.gml`, `.cnd`, `std` or `.prj` files in the Data Explorer or file conversion tools (e.g. `OGSFileConverter`):

<i class="far fa-exclamation-triangle"></i> Error message:

```bash
Error XSDError in http://www.opengeosys.org/images/xsd/OpenGeoSysCND.xsd, at line 1, column 1: Start tag expected.
Error XSDError in file:///../bc/well.cnd, at line 5, column 195: Loaded schema file is invalid.
XMLInterface::isValid() - XML File is invalid (in reference to schema ./OpenGeoSysCND.xsd).
Error XSDError in http://www.opengeosys.org/images/xsd/OpenGeoSysCND.xsd, at line 1, column 1: Start tag expected.
```

<i class="far fa-arrow-right"></i> Solution:

Open the affected file (e.g. `well.cnd` in this case) in a text editor and remove the following parameter of the XML root element (the first element in `< >`-brackets in the element, e.g. `<OpenGeoSysCND>`):

```xml
xsi:noNamespaceSchemaLocation="http://www.opengeosys.org/images/xsd/OpenGeoSysXXX.xsd"
```

where `XXX` can be `CND`, `GLI`, `PRJ` or `STN` depending on the file type.

Now the save the modified file and try to load it again.

<details>
    <summary>Background info:</summary>
    The XSD files may be downloaded from a web location. We changed the protocol of our web site to `https://` but due to some weird behaviour of the Qt XML validation code it tries to download the file (even if it is available locally) and does not respect the URL redirection to `https://` of the web server. Simply removing the part solves the problem. The XML is still validated! Newer OGS versions do not write that parameter into files anymore, see https://github.com/ufz/ogs/pull/2198.
</details>
