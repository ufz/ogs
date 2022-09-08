import tempfile
import os

import pytest

import OpenGeoSys


def test_module():
    arguments = [
        "",
        f"{os.path.abspath(os.path.dirname(__file__))}/../Data/Parabolic/LiquidFlow/Flux/cube_1e3_calculatesurfaceflux.prj",
        "-o " + tempfile.mkdtemp(),
    ]

    print("Python OpenGeoSys.init ...")
    OpenGeoSys.initialize(arguments)
    print("Python OpenGeoSys.executeSimulation ...")
    OpenGeoSys.executeSimulation()
    print("Python OpenGeoSys.finalize() ...")
    OpenGeoSys.finalize()


from . import push_argv


def _run(program, args):
    func = getattr(OpenGeoSys, program)
    args = ["%s.py" % program] + args
    with push_argv(args), pytest.raises(SystemExit) as excinfo:
        func()
    assert 0 == excinfo.value.code


def test_binaries():
    ignore_list = ["moveMeshNodes", "mpmetis", "tetgen"]  # have no --version cli flag
    for f in OpenGeoSys.binaries_list:
        if f not in ignore_list:
            _run(f, ["--version"])
