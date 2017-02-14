+++
date = "2017-01-14T22:56:13+01:00"
title = "Set Up Prerequisites"
author = "Lars Bilke"
weight = 2

[menu]
  [menu.devguide]
    parent = "getting-started"
    weight = 1
+++

## Minimum requirements

The minimum prerequisites to build OGS are:

- Git (version control tool, at least version 1.7.x)
- CMake (build configuration tool, at least version 3.1)
- A compiler with C++11-support
- Conan package manager (TODO crossref) OR required libraries (TODO crossref)

## Step: Install a compiler

### Windows

Because the C++11-standard is supported since **Visual C++ 2013** you will need at least **Windows 7** (64-bit recommended). It is perfectly fine to use the free Community Editions of Visual C++. We support the **2013 edition** and newer.

Download Visual Studio and install it:

- [Visual Studio Community 2015](https://go.microsoft.com/fwlink/?LinkId=532606&clcid=0x409)

### Linux

If you have a recent linux distribution you should also have a recent gcc. Please check that you have at least **gcc 4.9**:

{{< highlight bash >}}
$ gcc --version
gcc (GCC) 4.6.1
{{< /highlight >}}

### Mac

Please install Xcode from the App Store. Then please run the following command in the terminal to install the command line tools:

{{< highlight bash >}}
$ xcode-select --install
{{< /highlight >}}

Open Xcode one time to install some other Xcode stuff.

Now also install the Homebrew package manager:

{{< highlight bash >}}
$ ruby -e "$(curl -fsSL https://raw.github.com/mxcl/homebrew/go/install)"
$ brew doctor
{{< /highlight >}}

The Homebrew package manager is needed for installing other libraries and packages. It is just like a linux package manager.
