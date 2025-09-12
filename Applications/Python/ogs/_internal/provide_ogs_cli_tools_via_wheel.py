import os
import platform
import subprocess
import sys

from . import OGS_USE_PATH
from .binaries_list import binaries_list
from .get_bin_dir import get_bin_dir

OGS_BIN_DIR = get_bin_dir()
if platform.system() == "Windows":
    os.add_dll_directory(OGS_BIN_DIR)


# Not used when OGS_USE_PATH is true!
def ogs():
    raise SystemExit(ogs_with_args(sys.argv))


def ogs_with_args(argv):
    from ogs.OGSSimulator import (  # noqa: PLC0415
        OGSSimulation,
        check_command_line_arguments,
    )

    return_code = check_command_line_arguments(argv)
    # map mangled TCLAP status to usual exit status
    if return_code == 3:  # EXIT_ARGPARSE_FAILURE
        return 1  # EXIT_FAILURE
    if return_code == 2:  # EXIT_ARGPARSE_EXIT_OK
        return 0  # EXIT_SUCCESS

    if return_code != 0:
        return return_code

    try:
        sim = OGSSimulation(argv)
    except RuntimeError:
        print("Could not construct OGSSimulation object")
        return_code = 1  # EXIT_FAILURE
    else:
        return_code = sim.execute_simulation()
        sim.close()
    return return_code


def _program(name, args):
    exe = OGS_BIN_DIR / name
    env = None  # by default use unmodified environment
    if OGS_USE_PATH:
        exe = name
        env = os.environ.copy()
        # prevent infinite recursion if OGS in PATH happens to be this very
        # script
        env["OGS_USE_PATH"] = "0"
        print(f"OGS_USE_PATH is true: {name} from $PATH is used!")
    return subprocess.run([exe] + args, env=env).returncode  # noqa: PLW1510


FUNC_TEMPLATE = """def {0}(): raise SystemExit(_program("{0}", sys.argv[1:]))"""
for f in binaries_list:
    if f == "ogs" and not OGS_USE_PATH:
        continue  # provided by separate function
    # When OGS_USE_PATH is true then ogs()-function above is not used!
    exec(FUNC_TEMPLATE.format(f))
