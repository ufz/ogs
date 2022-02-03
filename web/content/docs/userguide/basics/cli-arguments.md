+++
date = "2021-05-06T13:00:13+01`:00"
title = "Command-line arguments"
author = "Lars Bilke"
weight = 3

[menu]
  [menu.userguide]
    parent = "basics"
+++

The following arguments are available:

```
$ ogs --help

   --enable-fpe
     enables floating point exceptions

   --unbuffered-std-out
     use unbuffered standard output

   --config-warnings-nonfatal
     warnings from parsing the configuration file will not trigger program
     abortion

   -l <LOG_LEVEL>,  --log-level <LOG_LEVEL>
     the verbosity of logging messages: none, error, warn, info, debug, all

   -o <PATH>,  --output-directory <PATH>
     the output directory to write to

   --write-prj
     Writes processed project file to output path /
     [prj_base_name]_processed.prj.

   -p <>,  --xml-patch <>  (accepted multiple times)
     the xml patch file(s) which is (are) applied (in the given order) to
     the PROJECT_FILE

   -r <PATH>,  --reference <PATH>
     Run output result comparison after successful simulation comparing to
     all files in the given path. This requires test definitions to be
     present in the project file.

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <PROJECT_FILE>
     (required)  Path to the ogs6 project file.
```
