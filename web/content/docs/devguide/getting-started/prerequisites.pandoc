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
- A compiler with [C++17](http://en.wikipedia.org/wiki/C%2B%2B17)-support
- [Conan package manager](https://www.conan.io/) (at least version {{< dataFile "versions.minimum_version.conan" >}}) **OR** install [required libraries]({{< ref "third-party-libraries.pandoc" >}}) manually (for advanced users only!)

## Step: Install a compiler

::: {.win}

::: {.note}

### Alternative setup

Please note that the following setup on Windows is the **native Windows development setup**. This native setup is **quite involved** and **heavy on system resources**. We can recommend an alternative setup in which the [Windows Subsystem for Linux](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux) is used: Setup and development of OGS follows the Linux way but you can use your Windows IDE (especially Visual Studio Code) for development and debugging. If this sounds interesting please [follow the steps here]({{< ref "wsl.pandoc" >}})!
:::

As we use lots of features of the C++17-standard we support **Visual Studio {{< dataFile "versions.minimum_version.msvc.year" >}}** and up. Therefore you will need at least **Windows 7** (64-bit recommended). It is perfectly fine to use the free Community Edition of Visual Studio.

- Download and install [Visual Studio Community](https://www.visualstudio.com)
  - Select the *workload* `Desktop Development with C++`
  - You can uncheck everything else
- When installation finished please start Visual Studio once (when asked for credentials enter your Microsoft account or click on **Skip for now**)
:::

::: {.linux}
On Debian-based (e.g. Ubuntu) you need to install the `build-essential`-package (which contains the `gcc`-compiler and the `make`-tool):

```bash
sudo apt install build-essential
```

You need to have at least **gcc {{< dataFile "versions.minimum_version.gcc" >}}**:

```bash
$ gcc --version
gcc (GCC) {{< dataFile "versions.minimum_version.gcc" >}}.0
```

::: {.note}

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

:::

:::

::: {.mac}
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
:::

## Step: Install Git

Git is a powerful and distributed version control system. OGS source code is hosted on [GitLab](https://gitlab.opengeosys.org/ogs/ogs). See the developer guide page on [Code Reviews]({{< ref "code-reviews" >}}) for more info on how OGS uses GitLab for collaborative development.

::: {.win}
Download and install git from the [git homepage](http://git-scm.com/download/win). Use the default installer options but also enable `Enable symbolic links` under the *Configuring extra options* page.

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

:::

::: {.linux}
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

:::

::: {.mac}
Install Git with Homebrew:

```bash
brew install git
```

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

:::

## Step: Install Python 3

::: {.win}

- [Download the Python 3 installer](https://www.python.org/ftp/python/3.7.2/python-3.7.2-amd64-webinstall.exe)
- Install with the following options
  - Check *Add Python 3.X to PATH*
  - *Customize installation*
  - Make sure to have `pip` enabled (you may uncheck *Documentation*, *tcl/tk*, *Python test suite*)
  - You may check *Install for all users*
  - Check *Add Python to environment variables*!

:::

::: {.linux}

Install Python 3 and pip:

```bash
sudo apt-get install python3 python3-pip
```

:::

::: {.mac}
Install Python 3 with Homebrew:

```bash
brew install python
```

:::

## Step: Install CMake

::: {.win}

- Download the installer, at the [CMake download page](http://www.cmake.org/cmake/resources/software.html) choose the **Windows (Win32 Installer)**.
- Execute the installer, please check the **Add CMake to the system path for all users**-option
:::

::: {.linux}
Install CMake via Kitware's APT Repository by [following their instructions](https://apt.kitware.com/).

For other linux distributions you want to use your distributions package manager, [pip](https://pypi.org/project/cmake/) or [snap](https://snapcraft.io/cmake).
:::

::: {.mac}
Install CMake with Homebrew:

```bash
brew install cmake
```

:::

## Step: Install Conan package manager

The [Conan package manager](https://www.conan.io) helps to install all required libraries in a convenient way on every platform. If you prefer you can also [install libraries manually]({{< ref "third-party-libraries.pandoc" >}}) instead.

Install Conan (>= {{< dataFile "versions.minimum_version.conan" >}}) with Python's pip:

```bash
pip3 install --user conan
```

This installed `conan` to `.local/bin` (or `C:\Users\[username]\AppData\Roaming\Python\Python37\Scripts` on Windows) in your home directory. Make sure to have this directory in your `PATH`!

Check on a newly opened command line if Conan was installed successfully:

```bash
$ conan --version
Conan version {{< dataFile "versions.minimum_version.conan" >}}
```
