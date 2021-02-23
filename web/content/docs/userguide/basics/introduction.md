+++
date = "2018-02-27T11:00:13+01:00"
title = "Introduction"
author = "Lars Bilke"
weight = 1

aliases = [ "/docs/userguide/",
            "/docs/quickstart/",
            "/docs/quickstart/basics/quickstart" ] # Redirect for Hydrology III tutorial

[menu.docs]
name = "User Guide"
identifier = "userguide"
weight = 1
post = "Download, install and run an OGS benchmark in 5 minutes! No development setup required."

[menu]
  [menu.userguide]
    parent = "basics"
+++

## Download

Download the latest release of OpenGeoSys from the [Releases](/releases)-page. Be sure to pick the correct file for your operating system.

## Installation

OGS itself is a simple executable file so you can put it anywhere you like. For convenience you may put into a location which is in your `PATH`-environment variable which allows you to start the executable without specifying its full file path.

## Download benchmarks

You can download the latest benchmark files from GitLab:

- On our OGS repository on GitLab browse to the [Tests/Data](https://gitlab.opengeosys.org/ogs/ogs/-/tree/master/Tests/Data)-folder
- Browse to the process subfolder you are interested in, e.g. [Elliptic](https://gitlab.opengeosys.org/ogs/ogs/-/tree/master/Tests/Data/Elliptic) (1.)
- Find the Downloads-dropdown (2.)
- Choose an appropriate format under **Download this directory** (3.)
- Uncompress the downloaded file

See [the Benchmarks section](/docs/benchmarks/) for more information on the benchmarks.

![Instructions for downloading benchmarks](/docs/userguide/basics/Download_Benchmarks.png)

## Running

OGS is a command line application and requires the path to a `.prj`-file as an argument.

<div class='win'>

To run it open a new command line shell (called *cmd.exe*). Now simply type `ogs` (if the executable is in your `PATH`-environment variable) or specify its full path (e.g.: `C:\Users\MyUserName\Downloads\ogs.exe`) and hit `ENTER`.

OGS prints out its usage instructions:

```bash
PARSE ERROR:
             Required argument missing: project-file

Brief USAGE:
   ogs  [--] [--version] [-h] <PROJECT FILE>

For complete USAGE and HELP type:
   ogs --help
```

You can see that there is the project-file missing.

Then simply supply the path to a project file as an argument to the OGS executable:

```bash
ogs .\Path\to\BenchmarkName.prj
```

</div>

<div class='linux'>

To run it open a new command line shell (*Terminal*). Now simply type `ogs` (if the executable is in your `PATH`-environment variable) or specify its full path (e.g.: `./path/to/ogs`) and hit `ENTER`.

OGS prints out its usage instructions:

```bash
PARSE ERROR:
             Required argument missing: project-file

Brief USAGE:
   ogs  [--] [--version] [-h] <PROJECT FILE>

For complete USAGE and HELP type:
   ogs --help
```

You can see that there is the project-file missing.

Then simply supply the path to a project file as an argument to the OGS executable:

```bash
ogs ./path/to/BenchmarkName.prj
```

</div>

<div class='mac'>

See Linux tab!

</div>
