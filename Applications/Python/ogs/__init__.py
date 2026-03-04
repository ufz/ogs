import importlib.util
import os
import sys
from pathlib import Path

# on wheel build config.py is not generated at the beginning of the build
# but later after CMake has run. Also cli works later when already build only.
config_spec = importlib.util.find_spec(".config", package="ogs")
if config_spec is not None:
    from ._internal.wrap_cli_tools import cli  # noqa: F401
    from .config import *  # noqa: F403

# DLL loading in Python does not respect the PATH, see
# https://docs.python.org/3/library/os.html#os.add_dll_directory
if sys.platform == "win32" and (root := os.environ.get("MKLROOT")):
    os.add_dll_directory(str(Path(root) / "bin"))
