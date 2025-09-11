import importlib.util

# on wheel build config.py is not generated at the beginning of the build
# but later after CMake has run. Also cli works later when already build only.
config_spec = importlib.util.find_spec(".config", package="ogs")
if config_spec is not None:
    from ._internal.wrap_cli_tools import cli  # noqa: F401
    from .config import *  # noqa: F403
