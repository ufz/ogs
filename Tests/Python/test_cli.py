import tempfile
import os

import pytest

import ogs

from . import push_argv


def _run(program, args):
    func = getattr(ogs, program)
    args = ["%s.py" % program] + args
    with push_argv(args), pytest.raises(SystemExit) as excinfo:
        func()
    assert 0 == excinfo.value.code


def test_binaries():
    ignore_list = ["moveMeshNodes", "mpmetis", "tetgen"]  # have no --version cli flag
    for f in ogs.binaries_list:
        if f not in ignore_list:
            _run(f, ["--version"])
