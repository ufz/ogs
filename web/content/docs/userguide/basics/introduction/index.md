+++
date = "2018-02-27T11:00:13+01:00"
title = "Introduction"
author = "Lars Bilke"
weight = 12
+++

## Installation

There are various ways to obtain a (pre-built) running version of OpenGeoSys (OGS) on your machine. You can have OGS for different
operating systems including Windows, Linux, and macOS. The basic modelling platform is available for all operating systems.
The different operating systems and installation methods give you the feature matrix:

| Operating system / Installation method                                                                                       | [Processes](/docs/userguide/blocks/processes/)                                    | [MFront](/docs/userguide/features/mfront/) | [PETSc]({{< ref "parallel_computing_mpi" >}})
| ---------------------------------------------------------------------------------------------------------------------------- | :-------------------------------------------------------------------------------: | :----------------------------------------: | :-------------------------------------------:
| <i class="fab fa-windows"></i> Windows / [pip](#install-via-pip) & [binary download](#alternative-install-via-binary-downloads) | All |                     ❌                      |                       ❌
| <i class="fab fa-apple"></i> macOS / [pip](#install-via-pip)                                                                    |                                        All                                        |                     ✅                      |                       ❌
| <i class="fab fa-linux"></i> Linux / [pip](#install-via-pip) & [Serial Container]({{< ref "container.md" >}})                   |                                        All                                        |                     ✅                      |                       ❌
| <i class="fab fa-linux"></i> Linux / [PETSc Container]({{< ref "container.md" >}})                                              |                                        All                                        |                     ✅                      |                       ✅

<div class="note">

### Using Linux binaries on other operating systems

Please note that you can use Linux binaries installed via `pip` or in the form of a container also on other operating systems.

For Windows we recommend the [Windows Subsystem for Linux (WSL)]({{< ref "wsl" >}}).

On macOS you can use either a virtual machine (e.g. via [UTM](https://docs.getutm.app/installation/macos/)) or run a [Docker container]({{< ref "container.md#with-docker" >}}).

</div>

## Install via pip

A straightforward way of installing OGS is via Python's [`pip`-tool](https://packaging.python.org/en/latest/tutorials/installing-packages/):

<div class='win'>

```powershell
# Optional: create a Python virtual environment, see below
python -m venv .venv
.\venv\Scripts\Activate.ps1 # for PowerShell OR
.\venv\Scripts\activate.bat # for Command Prompt


# Install latest ogstools and ogs release ({{< ogs-last-release >}}) via pip:
pip install ogstools[ogs]

# OR: Install latest ogs master (weekly) via pip:
pip install --pre --index-url https://gitlab.opengeosys.org/api/v4/projects/120/packages/pypi/simple ogstools[ogs]
```

If you get errors when calling the `Activate.ps1` script you have to run the following command once:

```powershell
Set-Execution-Policy Unrestricted -Scope CurrentUser
```

</div>

<div class='linux'>

```bash
# Optional: create a Python virtual environment, see below
python -m venv .venv      # or `python3 -m venv .venv`
source .venv/bin/activate

# Install latest ogstools and ogs release ({{< ogs-last-release >}}) via pip:
pip install ogstools[ogs]

# OR: Install latest ogs master (weekly) via pip:
pip install --pre --index-url https://gitlab.opengeosys.org/api/v4/projects/120/packages/pypi/simple ogstools[ogs]
```

</div>

<div class='mac'>
See Linux-tab!
</div>

Please note that this requires a Python installation (version 3.10 - 3.13) with the `pip`-tool.

We recommend using Python within a [virtual environment](https://docs.python.org/3/library/venv.html) to keep possible
conflicts of different Python-packages localised. If you use `pip` for installation of OGS in a virtual environment and you
activate the virtual environment, then OGS and its tools are automatically also in the `PATH`. If the virtual environment is
not activated you may still use OGS, but either have to give the full path to `ogs` being located in the `bin` folder (or `Scripts` folder on Windows) of the
virtual environment, or add this path to your `PATH`-environment. Moreover, `pip` may print instructions which directory needs
 to be added to the `PATH`.

<div class="note">

### Advanced usage: Use binaries from `PATH` instead of pip package

When using the OGS Python bindings it may be desired to use a different OGS version than provided by the pip package. This can be achieved by setting the environment variable `OGS_USE_PATH`. E.g.:

```bash
pip install ogs
# now ogs from pip is used
export PATH=/path/to/ogs/bin:$PATH
OGS_USE_PATH=1 python
# now in Python interpreter:
>>> import ogs
>>> ogs.cli.ogs("--version") # uses /path/to/ogs/bin/ogs instead of ogs from pip
```

</div>

<div class='win'>

<div class="note">

### Alternative: Install via binary downloads

Another way to obtain a running version is
to just download the latest stable or development release of OpenGeoSys from the [Releases](/releases)-page.

By downloading from the release page you will get a bunch of folders and files. However, OGS itself will come as a simple
executable file, which you will find in the `bin` sub-folder. You can put the executable wherever you like. For convenience you
may put it into a location which is in your `PATH`-environment variable or [modify](https://stackoverflow.com/a/714918/80480) the `PATH` environment variable which allows you to start the executable without
specifying its full file path. Then just calling `ogs` from the terminal is sufficient.

</div>

</div>

<div class='linux'>

</div>

<div class='mac'>

</div>

## Download benchmarks

You can download the latest benchmark files from GitLab:

- On our OGS repository on GitLab browse to the [Tests/Data](https://gitlab.opengeosys.org/ogs/ogs/-/tree/master/Tests/Data)
-folder
- Browse to the process subfolder you are interested in, e.g. [Elliptic](https://gitlab.opengeosys.org/ogs/ogs/-/tree/master/Tests/Data/Elliptic) (1.)
- Find the Downloads-dropdown (2.)
- Choose an appropriate format under **Download this directory** (3.)
- Uncompress the downloaded file

See [the Benchmarks section](/docs/benchmarks/) for more information on the benchmarks.

![Instructions for downloading benchmarks](/docs/userguide/basics/Download_Benchmarks.png)

## Running

OGS is a command line application and requires the path to a `.prj`-file as an argument.

<div class='win'>

To run it open a new command line shell (called *cmd.exe*). Now simply type `ogs` (if the executable is in your `PATH`
-environment variable) or specify its full path (e.g.: `C:\Users\MyUserName\Downloads\ogs.exe`) and hit `ENTER`.

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

To run it open a new command line shell (*Terminal*). Now simply type `ogs` (if the executable is in your `PATH`-environment
variable) or specify its full path (e.g.: `./path/to/ogs`) and hit `ENTER`.

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

To run it open a new command line shell (*Terminal*). Now simply type `ogs` (if the executable is in your `PATH`-environment
variable) or specify its full path (e.g.: `./path/to/ogs`) and hit `ENTER`.

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
