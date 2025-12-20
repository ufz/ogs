# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

import ogs._internal.provide_ogs_cli_tools_via_wheel as ogs_cli_wheel
import pytest
from ogs._internal.binaries_list import binaries_list

from . import push_argv


def _run(program, args):
    func = getattr(ogs_cli_wheel, program)
    args = [f"{program}.py"] + args
    with push_argv(args), pytest.raises(SystemExit) as excinfo:
        func()
    assert excinfo.value.code == 0


def test_binaries():
    ignore_list = ["moveMeshNodes", "mpmetis"]  # have no --version cli flag
    for f in binaries_list:
        if f not in ignore_list:
            _run(f, ["--version"])
