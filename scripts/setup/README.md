# Setup scripts for Windows #

## setup.sh ##

If you are on Windows you can run `setup.sh` on the *Git Bash* to setup all required third party libraries for OGS (including the Data Explorer). Make also sure you have working Internet connection. For more details have a look on the [devguide][dev].

The following options are available:

- `-a [x32|x64]` to build for either 32- or 64-bits
- `-s` to build Qt with Oracle SQL bindings

A typical call would look like this:

		cd /scripts/setup
		./setup.sh -a x64 -d

The compilation of the libraries will take some time. Make sure to set the Qt environment variables as described in the [devguide][dev].

[dev]: (http://ufz.github.com/devguide/win-prerequisites/).