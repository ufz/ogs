import importlib.util
import os

from ._internal.wrap_cli_tools import cli  # noqa: F401

if "PEP517_BUILD_BACKEND" in os.environ:
    # on wheel build config.py is not generated at the beginning of the build
    # but later after CMake has run
    config_spec = importlib.util.find_spec(".config", package="ogs")
    if config_spec is not None:
        from .config import *  # noqa: F403
else:
    from .config import *  # noqa: F403
