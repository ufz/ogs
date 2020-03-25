+++
date = "2018-02-26T11:00:13+01:00"
title = "Log and Debug Output"
author = "Lars Bilke"
weight = 1034

[menu]
  [menu.devguide]
    parent = "advanced"
+++

## Introduction

For application output we use [spdlog](https://github.com/gabime/spdlog) which
is already integrated in OGS. Spdlog provides several verbosity levels which can
be used with simple calls:

```cpp
ERR("An error message!");
WARN("A warning message.");
INFO("An information message...");
DBUG("A debug message.");
```

As arguments you can use the same functionality as in [fmt](https://fmt.dev)---a
modern formatting library:

```cpp
int foo = 42;
double boo = 3.14;
WARN("Foo is {}! Current value is {:10.2g}.", foo, boo);
```

For more information see the [spdlog
wiki](https://github.com/gabime/spdlog/wiki).

On release builds the default log level is `INFO`, for debug builds it is
`DEBUG`.
The log level can be adjusted on the command line with `-l <LOG_LEVEL>` flag.
