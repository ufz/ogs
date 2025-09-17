import os
import sysconfig
from pathlib import Path


def get_bin_dir():
    site_packages_path = sysconfig.get_paths()["purelib"]
    OGS_BIN_DIR = Path(site_packages_path) / "bin"

    if not OGS_BIN_DIR.exists():
        # bin in regular cmake build-directory
        if "CI" not in os.environ:
            print(
                "Warning: OGS_BIN_DIR does not exist, falling back to possible "
                "build directory location."
            )
        OGS_BIN_DIR = Path(__file__).parent.parent.parent.parent / "bin"

    return OGS_BIN_DIR
