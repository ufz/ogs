import importlib.util
import json
import sysconfig
from importlib.metadata import Distribution, PackageNotFoundError
from pathlib import Path


def get_bin_dir():
    pkg_is_editable = False
    if importlib.util.find_spec("ogs") is not None:
        try:
            dist = Distribution.from_name("ogs")
            direct_url = dist.read_text("direct_url.json")
            if direct_url:
                pkg_is_editable = (
                    json.loads(direct_url).get("dir_info", {}).get("editable", False)
                )
        except PackageNotFoundError:
            # ogs is importable but not an installed distribution (e.g. source tree, dev mode)
            pkg_is_editable = False
        except FileNotFoundError:
            # distribution is found but no direct_url.json
            pass

    if pkg_is_editable:
        site_packages_path = sysconfig.get_paths()["purelib"]
        OGS_BIN_DIR = Path(site_packages_path) / "bin"
    else:
        # Here, we assume that this script is installed, e.g., in a virtual environment
        # alongside a "bin" directory.
        OGS_BIN_DIR = Path(__file__).parent.parent.parent / "bin"  # installed wheel

    if not OGS_BIN_DIR.exists():
        # bin in regular cmake build-directory
        print(
            "Warning: OGS_BIN_DIR does not exist, falling back to possible "
            "build directory location."
        )
        OGS_BIN_DIR = Path(__file__).parent.parent.parent.parent / "bin"

    return OGS_BIN_DIR
