+++
date = "2017-01-14T22:56:13+01:00"
title = "Set Up Prerequisites"
author = "Lars Bilke"
weight = 2

[menu]
  [menu.devguide]
    parent = "getting-started"
+++

## Minimum requirements

The minimum prerequisites to build OGS are:

- An 64-bit operating system (Linux, Windows 7 and up, macOS)
- Git (version control tool, at least version {{< dataFile "versions.minimum_version.git" >}})
- CMake (build configuration tool, at least version {{< dataFile "versions.minimum_version.cmake" >}})
- A compiler with [C++20](http://en.wikipedia.org/wiki/C%2B%2B20)-support
- Python interpreter and libraries (and optionally the `uv` virtual environment management tool)
- *Optional (but recommended)*: [Ninja](https://ninja-build.org) build tool

<div class='note'>

### Note about skipping installation steps

A fresh system with none of the prerequisites fulfilled is assumed. Skipping installation steps or using a non-supported version might result in unexpected problems. If possible, you may consider reinstalling or manual modifying the configuration of the already installed tool.

</div>

## Step: Install a compiler

<div class='win'>

<div class='note'>

### Alternative setup

Please note that the following setup on Windows is the **native Windows development setup**. This native setup is **quite involved** and **heavy on system resources**. We can recommend an alternative setup in which the [Windows Subsystem for Linux](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux) is used: Setup and development of OGS follows the Linux way but you can use your Windows IDE (especially Visual Studio Code) for development and debugging. If this sounds interesting please [follow the steps here]({{< ref "wsl.md" >}})!

</div>

As we use lots of features of the C++20-standard we support **Visual Studio {{< dataFile "versions.minimum_version.msvc.year" >}}** with compiler version **{{< dataFile "versions.minimum_version.msvc.compiler" >}}** and up. Therefore you will need at least **Windows 7** (64-bit required). It is perfectly fine to use the free Community Edition of Visual Studio.

- Download and install [Visual Studio Community](https://www.visualstudio.com)
  - Select the *workload* `Desktop Development with C++`
  - You can uncheck everything else
- When installation finished please start Visual Studio once (when asked for credentials enter your Microsoft account or click on **Skip for now**)

</div>

<div class='linux'>

On Debian-based (we recommend using Ubuntu {{< dataFile "versions.tested_version.ubuntu" >}}) you need to install the `build-essential`-package (which contains the `gcc`-compiler and the `make`-tool):

```bash
sudo apt install build-essential
```

You need to have at least **GCC {{< dataFile "versions.minimum_version.gcc" >}}** which you can check with `gcc --version`.

<div class='note'>

### Install the required compiler on older Ubuntu versions

If you are on an older Ubuntu version than {{< dataFile "versions.tested_version.ubuntu" >}} you can install a newer compiler from the `ubuntu-toolchain-r/test`-repository (with the following steps e.g. you can install GCC {{< dataFile "versions.minimum_version.gcc" 1 >}} on Ubuntu 20.04):

```bash
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt-get install gcc-{{< dataFile "versions.minimum_version.gcc" 1 >}}
sudo apt-get install g++-{{< dataFile "versions.minimum_version.gcc" 1 >}}
```

To make the newly installed compiler the default one:

```bash
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-{{< dataFile "versions.minimum_version.gcc" 1 >}} 60 \
  --slave /usr/bin/g++ g++ /usr/bin/g++-{{< dataFile "versions.minimum_version.gcc" 1 >}}
```

If you do not do this you have to specify the compiler during the first CMake run:

```bash
CC=gcc-{{< dataFile "versions.minimum_version.gcc" 1 >}} CXX=c++-{{< dataFile "versions.minimum_version.gcc" 1 >}} cmake ../ogs [more CMake options]
```

</div>

</div>

<div class='mac'>
Please install Xcode from the App Store. Then please run the following command in the terminal to install the command line tools:

```bash
xcode-select --install
```

Open Xcode one time to install some other Xcode stuff.

Now also install the [Homebrew](http://brew.sh) package manager:

```bash
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew doctor
```

The Homebrew package manager is needed for installing other libraries and packages. It is just like a Linux package manager.
</div>

## Step: Install Git

Git is a powerful and distributed version control system. OGS source code is hosted on [GitLab](https://gitlab.opengeosys.org/ogs/ogs). See the developer guide page on [Code Reviews]({{< ref "code-reviews" >}}) for more info on how OGS uses GitLab for collaborative development.

<div class='win'>

**Download** and install git from the [git homepage](http://git-scm.com/download/win). Use the default installer options but also enable `Enable symbolic links` under the *Configuring extra options* page.

![Enable symbolic links option](git-installer-win.png)

This install a new command line called *Git Bash* which should be used for all git operations.

Let Git know who you are:

```bash
git config --global user.name "Your Name Here"
git config --global user.email "your_email@example.com"
```

Enable the long paths feature:

```bash
git config --global core.longpaths true
```

You may also have [enable the long paths feature in the Windows registry](https://www.thewindowsclub.com/how-to-enable-or-disable-win32-long-paths-in-windows-11-10).

In some corporate environments you may have to use a proxy server. In this case tell git about it:

```bash
git config --global http.proxy http://yourproxy.example.com
```

</div>

<div class='linux'>

Please check if Git is already installed:

```bash
$ git --version
git version {{< dataFile "versions.minimum_version.git" >}}
```

Otherwise please install Git with your favorite package manager:

```bash
sudo apt-get install git
```

Let Git know who you are:

```bash
git config --global user.name "Your Name Here"
git config --global user.email "your_email@example.com"
```

Optionally enable password storing when interacting with a remote server:

```bash
git config --global credential.helper store
```

In some corporate environments you may have to use a proxy server. In this case tell git about it:

```bash
git config --global http.proxy http://yourproxy.example.com
```

</div>

<div class='mac'>

Git is already installed.

Let Git know who you are:

```bash
git config --global user.name "Your Name Here"
git config --global user.email "your_email@example.com"
```

The [graphical GitHub client](http://mac.github.com/) is also maybe worth a look.

In some corporate environments you may have to use a proxy server. In this case tell git about it:

```bash
git config --global http.proxy http://yourproxy.example.com
```

</div>

## Step: Install CMake

<div class='win'>

- Download the installer, at the [CMake download page](http://www.cmake.org/cmake/resources/software.html) choose the **Windows (Win32 Installer)**.
- Execute the installer, please check the **Add CMake to the system path for all users**-option

</div>

<div class='linux'>

Install CMake using system package manager, *e.g.* on Ubuntu

```bash
sudo apt install cmake
```

For other Linux distributions you want to use your distributions package manager, [pip](https://pypi.org/project/cmake/) or [snap](https://snapcraft.io/cmake).

</div>

<div class='mac'>

Install CMake with Homebrew:

```bash
brew install cmake
```

</div>

## Optional: Install Ninja

We recommend [`ninja`](https://ninja-build.org) as a cross-platform build tool (`make`-replacement).

<div class='win'>

Download the [binary from GitHub](https://github.com/ninja-build/ninja/releases) and put the extracted `ninja.exe` in the `PATH`.

</div>

<div class='linux'>

```bash
sudo apt-get install ninja-build
```

</div>

<div class='mac'>

Install Ninja with Homebrew:

```bash
brew install ninja
```

</div>

## Install Python 3

<div class='win'>

- [Download the Python 3 installer](https://www.python.org/ftp/python/3.10.10/python-3.10.10-amd64.exe)
- Install with the following options
  - Check *Add Python 3.X to PATH*
  - *Customize installation*
  - Make sure to have `pip` enabled (you may uncheck *Documentation*, *tcl/tk*, *Python test suite*)
  - Check *Add Python to environment variables*!

</div>

<div class='linux'>

Install Python 3 and pip:

```bash
sudo apt-get install python3 python3-pip
```

</div>

<div class='mac'>

Install Python 3 with Homebrew:

```bash
brew install python
```

</div>

If you want to execute or develop Jupyter notebooks for more interactive benchmarks please install [`uv`](https://docs.astral.sh/uv/getting-started/installation/), a tool which handles Python dependencies and virtual environments. See the [Jupyter page]({{< ref "jupyter-docs.md" >}}) for more info.

## Optional: Install Qt, NetCDF and other dependencies for the Data Explorer

Use [Another Qt installer (`aqt`)](https://github.com/miurahr/aqtinstall) for installing the Qt binaries to some path on your machine:

<div class='win'>

```bat
pip install aqtinstall
mkdir qt
cd qt
aqt install-qt win desktop {{< dataFile "versions.tested_version.qt" >}} win64_msvc2019_64
aqt install-qt win desktop {{< dataFile "versions.tested_version.qt" >}} win64_msvc2019_64 --archives qtxmlpatterns
```

This will install Qt to `[your-directory]/qt/{{< dataFile "versions.tested_version.qt" >}}/msvc2019_64`.

To finish add `[your-directory]/qt/{{< dataFile "versions.tested_version.qt" >}}/msvc2019_64/bin` bin to the `PATH` environment variable.

Install **NetCDF4** by downloading and installing the [official installer](https://downloads.unidata.ucar.edu/netcdf-c/4.9.2/netCDF4.9.2-NC4-64.exe). The C++-bindings to NetCDF are automatically build via CPM.

</div>

<div class='linux'>

```bash
pip install aqtinstall
mkdir /opt/qt
cd /opt/qt
aqt install-qt linux desktop {{< dataFile "versions.tested_version.qt" >}} gcc_64
aqt install-qt linux desktop {{< dataFile "versions.tested_version.qt" >}} gcc_64 --archives qtxmlpatterns qtx11extras
```

Make sure to add `/opt/qt/{{< dataFile "versions.tested_version.qt" >}}/gcc_64/bin` to the `PATH`.

Install more dependencies for VTK rendering and for NetCDF IO:

```bash
sudo apt-get install freeglut3 freeglut3-dev libglew-dev libglu1-mesa libglu1-mesa-dev \
  libgl1-mesa-glx libgl1-mesa-dev libnetcdf-c++4-dev
```

</div>

<div class='mac'>

```bash
pip install aqtinstall
mkdir /opt/qt
cd /opt/qt
aqt install-qt mac desktop {{< dataFile "versions.tested_version.qt" >}} clang_64
aqt install-qt mac desktop {{< dataFile "versions.tested_version.qt" >}} clang_64 --archives qtxmlpatterns qtx11extras
```

Make sure to add `/opt/qt/{{< dataFile "versions.tested_version.qt" >}}/clang_64/bin` to the `PATH`.

Install NetCDF:

```bash
brew install netcdf-cxx
```

</div>
