import os
import platform
import sys

from ._internal.wrap_cli_tools import cli  # noqa: F401

if platform.system() != "Windows":
    # Otherwise runtime undefined symbols, see
    # https://stackoverflow.com/a/60841073/80480
    sys.setdlopenflags(os.RTLD_GLOBAL | os.RTLD_LAZY)
