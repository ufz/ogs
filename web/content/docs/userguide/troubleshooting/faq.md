+++
date = "2018-02-27T11:00:13+01:00"
title = " Frequently asked questions"
author = "Lars Bilke and Feliks Kiszkurno"
weight = 1
+++

## Windows only: `ogs`` or any of its tools cannot open an existing file with a long path

Windows has a [default path length limit of 260 characters](https://learn.microsoft.com/en-us/windows/win32/fileio/maximum-file-path-limitation). Especially in workflows this limit can be exceeded easily. To enable long paths on Windows you need to alter the Windows registry which requires administrative user privileges:

- Open a PowerShell command prompt as an Administrator
- Run the following script:

  ```powershell
  New-ItemProperty -Path "HKLM:\SYSTEM\CurrentControlSet\Control\FileSystem" -Name "LongPathsEnabled" -Value 1 -PropertyType DWORD -Force
  ```

- You may have to restart the computer

## `XSDError: Loaded schema file is invalid` error encountered when running DataExplorer

You may encounter the following error (or similar) on opening `.gml`, `.cnd`, `std` or `.prj` files in the Data Explorer or file
conversion tools (e.g., `OGSFileConverter`):

<i class="far fa-exclamation-triangle"></i> Error message:

```bash
Error XSDError in http://www.opengeosys.org/stable/images/xsd/OpenGeoSysCND.xsd, at line 1, column 1: Start tag expected.
Error XSDError in file:///../bc/well.cnd, at line 5, column 195: Loaded schema file is invalid.
XMLInterface::isValid() - XML File is invalid (in reference to schema ./OpenGeoSysCND.xsd).
Error XSDError in http://www.opengeosys.org/stable/images/xsd/OpenGeoSysCND.xsd, at line 1, column 1: Start tag expected.
```

<i class="far fa-arrow-right"></i> Solution:

Open the affected file (e.g. `well.cnd` in this case) in a text editor and remove the following parameter of the XML root element (the first element in `< >`-brackets in the element, e.g. `<OpenGeoSysCND>`):

```xml
xsi:noNamespaceSchemaLocation="http://www.opengeosys.org/stable/images/xsd/OpenGeoSysXXX.xsd"
```

where `XXX` can be `CND`, `GLI`, `PRJ` or `STN` depending on the file type.

Now, save the modified file and try to load it again.

<!-- vale off -->

<details>
    <summary>Background info:</summary>
    The XSD files may be downloaded from a web location. We changed the protocol of our web site to `https://` but due to some weird behavior of the Qt XML validation code it tries to download the file (even if it is available locally) and does not respect the URL redirection to `https://` of the web server. Simply removing the part solves the problem. The XML is still validated! Newer OGS versions do not write that parameter into files anymore, see [!2198](https://github.com/ufz/ogs/pull/2198).
</details>

## Oh, no! It diverges

OpenGeoSys will sometimes crash with an error message stating, that it couldn't reduce the time step any more. There are several
things that can be done to solve this:

- In [Time loop](/docs/userguide/blocks/time_loop/) adjust time stepping allowing smaller time steps to be taken.
- In [Time loop](/docs/userguide/blocks/time_loop/) adjust the [error tolerance](/docs/userguide/blocks/time_loop/#error-tolerances)
- Check if [mesh quality](/docs/userguide/blocks/meshes/#mesh-quality) is sufficient
- Verify if all [parameters](/docs/userguide/blocks/parameters/), especially those that vary during the simulation, make physical sense throughout the whole experiment and combined with other parameters.

## Oh, no! There are NaNs everywhere

If in the output log all values are [NaN](/docs/userguide/troubleshooting/glossary/#nan) and one or more parameters were defined
as functions, there is a good chance that a division by zero, or other undefined operation happened.
