+++
date = "2017-01-14T22:56:13+01:00"
title = "Set Up Prerequisites"
author = "Lars Bilke"
weight = 1002

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
- *Optional (but recommended)*: [Ninja](https://ninja-build.org) build tool
- *Optional*: [Conan package manager](https://www.conan.io/) (at least version {{< dataFile "versions.minimum_version.conan" >}}) for some optional dependencies.

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

As we use lots of features of the C++17-standard we support **Visual Studio {{< dataFile "versions.minimum_version.msvc.year" >}}** with compiler version **{{< dataFile "versions.minimum_version.msvc.compiler" >}}** and up. Therefore you will need at least **Windows 7** (64-bit required). It is perfectly fine to use the free Community Edition of Visual Studio.

- Download and install [Visual Studio Community](https://www.visualstudio.com)
  - Select the *workload* `Desktop Development with C++`
  - You can uncheck everything else
- When installation finished please start Visual Studio once (when asked for credentials enter your Microsoft account or click on **Skip for now**)

</div>

<div class='linux'>
On Debian-based (e.g. Ubuntu) you need to install the `build-essential`-package (which contains the `gcc`-compiler and the `make`-tool):

```bash
sudo apt install build-essential
```

You need to have at least **gcc {{< dataFile "versions.minimum_version.gcc" >}}**:

```bash
$ gcc --version
gcc (GCC) {{< dataFile "versions.minimum_version.gcc" >}}.0
```

<div class='note'>

### Install a newer compiler on Ubuntu

We recommend using Ubuntu {{< dataFile "versions.tested_version.ubuntu" >}} as its standard `gcc` package is already at version 9. If you are on an older Ubuntu version you can install a newer compiler from the `ubuntu-toolchain-r/test`-repository:

```bash
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt-get install gcc-9 g++-9
```

To make the newly installed compiler the default one:

```bash
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 60 \
  --slave /usr/bin/g++ g++ /usr/bin/g++-9
```

If you do not do this you have to specify the compiler during the first CMake run:

```bash
CC=gcc-9 CXX=c++-9 cmake ../ogs [more CMake options]
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

The Homebrew package manager is needed for installing other libraries and packages. It is just like a linux package manager.
</div>

## Step: Install Git

Git is a powerful and distributed version control system. OGS source code is hosted on [GitLab](https://gitlab.opengeosys.org/ogs/ogs). See the developer guide page on [Code Reviews]({{< ref "code-reviews" >}}) for more info on how OGS uses GitLab for collaborative development.

<div class='win'>

**Download** and install git from the [git homepage](http://git-scm.com/download/win). Use the default installer options but also enable `Enable symbolic links` under the *Configuring extra options* page.

![Enable symbolic links option](../git-installer-win.png)

This install a new command line called *Git Bash* which should be used for all git operations.

Let Git know who you are:

```bash
git config --global user.name "Your Name Here"
git config --global user.email "your_email@example.com"
```

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

Install CMake via Kitware's APT Repository by [following their instructions](https://apt.kitware.com/).

For other linux distributions you want to use your distributions package manager, [pip](https://pypi.org/project/cmake/) or [snap](https://snapcraft.io/cmake).
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

## Optional: Install Python 3

<div class='win'>

- [Download the Python 3 installer](https://www.python.org/ftp/python/3.9.1/python-3.9.1-amd64.exe)
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

## Optional: Install Qt for the Data Explorer

Use [Another Qt installer(aqt)](https://github.com/miurahr/aqtinstall) for installing the Qt binaries to some path on your machine:

<div class='win'>

```bat
pip install aqtinstall
mkdir qt
cd qt
aqt install {{< dataFile "versions.tested_version.qt" >}} windows desktop win64_msvc2019_64 -m qtxmlpatterns
```

This will install Qt to `[your-directory]/qt/{{< dataFile "versions.tested_version.qt" >}}/msvc2019_64`.

To finish add `[your-directory]/qt/{{< dataFile "versions.tested_version.qt" >}}/msvc2019_64/bin` bin to the `PATH` environment variable.

</div>

<div class='linux'>

```bash
pip install aqtinstall
mkdir /opt/qt
cd /opt/qt
aqt install {{< dataFile "versions.tested_version.qt" >}} linux desktop gcc_64 -m xmlpatterns,x11extras
```

Make sure to add `/opt/qt/{{< dataFile "versions.tested_version.qt" >}}/gcc_64/bin` to the `PATH`.

</div>

<div class='mac'>

```bash
pip install aqtinstall
mkdir /opt/qt
cd /opt/qt
aqt install {{< dataFile "versions.tested_version.qt" >}} linux desktop clang_64 -m xmlpatterns,x11extras
```

Make sure to add `/opt/qt/{{< dataFile "versions.tested_version.qt" >}}/clang_64/bin` to the `PATH`.

</div>

## Optional: Install Conan package manager

You only need Conan if you intend to build with one of the following settings **and** do not want to install their dependencies manually:

- `OGS_USE_NETCDF` – NetCDF IO, requires netcdf-cxx4

Install Conan (>= {{< dataFile "versions.minimum_version.conan" >}}) with Python's pip:

```bash
pip3 install --user conan
```

This installed `conan` to `.local/bin` (or `C:\Users\[username]\AppData\Roaming\Python\Python39\Scripts` on Windows) in your home directory. Make sure to have this directory in your `PATH`!

Check on a newly opened command line if Conan was installed successfully:

```bash
$ conan --version
Conan version {{< dataFile "versions.minimum_version.conan" >}}
```

<div class='note'>

**Advanced usage:** You can also have Conan auto-installed when using the CMake-option `OGS_USE_CONAN=auto`. See the page on [Python environment]({{< ref "python-env.md" >}}) for details.

</div>
