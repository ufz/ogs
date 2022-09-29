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
+++

## Installation

<div class='win'>

Download the latest release of OpenGeoSys from the [Releases](/releases)-page. Be sure to pick the correct file for your operating system.

OGS itself is a simple executable file so you can put it anywhere you like. For convenience you may put into a location which is in your `PATH`-environment variable which allows you to start the executable without specifying its full file path.

<div class="note">

### Alternative: Install via `pip`

You can also install ogs via Python's [`pip`-tool](https://packaging.python.org/en/latest/tutorials/installing-packages/):

```bat
pip install ogs
```

If you install into an activated [virtual environment](https://docs.python.org/3/library/venv.html) then ogs and its tools are automatically also in the `PATH`. Otherwise `pip` will print instructions which directory needs to be added to the `PATH`.

</div>

</div>

<div class='linux'>

Install via Python's [`pip`-tool](https://packaging.python.org/en/latest/tutorials/installing-packages/):

```bash
pip install ogs
```

You may want to set up and activate a [virtual environment](https://docs.python.org/3/library/venv.html) before.

You could also use [`pipx`](https://pypa.github.io/pipx/) to install into an isolated environment.

</div>

<div class='mac'>

See Linux tab!

</div>

<div class="note">

### Limitations of the `pip`-based installation

- Serial configuration only! For PETSc-support please use a [Singularity container]({{< relref "container" >}}).
- No embedded Python interpreter, i.e. no Python boundary conditions!
- A Python (3.8 - 3.11) installation with `pip` is required.

</div>

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
