import importlib.util
import os

# on wheel build config.py is not generated at the beginning of the build
# but later after CMake has run. Also cli works later when already build only.
config_spec = importlib.util.find_spec(".config", package="ogs")
if config_spec is not None:
    from .. import config

    if config.OGS_GUIX_BUILD == "ON":
        OGS_USE_PATH = True
    else:
        OGS_USE_PATH = os.getenv("OGS_USE_PATH", "False").lower() in ("true", "1", "t")
