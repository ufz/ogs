import os
import subprocess
import sys

from ._cli import *

OGS_BIN_DIR = os.path.join(os.path.join(os.path.dirname(__file__), "bin"))


def _program(name, args):
    return subprocess.call([os.path.join(OGS_BIN_DIR, name)] + args)


def ogs():
    raise SystemExit(_program("ogs", sys.argv[1:]))
