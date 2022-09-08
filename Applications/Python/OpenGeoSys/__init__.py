import os
import subprocess
import sys

try:
    from .lib64._cli import *
except ImportError:
    try:
        from .lib._cli import *
    except ImportError:
        print("ERROR: could not import OpenGeoSys Python module!")

OGS_BIN_DIR = os.path.join(os.path.join(os.path.dirname(__file__), "bin"))


def _program(name, args):
    return subprocess.call([os.path.join(OGS_BIN_DIR, name)] + args)


# Binary entrypoints
def ogs():
    raise SystemExit(_program("ogs", sys.argv[1:]))
