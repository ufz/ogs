# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

import os
import signal
import subprocess

import pytest


# In the past using OGS_USE_PATH=1 together with OGS installed via pip could
# lead to an infinite recursion when calling ogs or OGS utilities.
# This test checks that no such infinite recursion occurs anymore.
#
# It's not a completely watertight test, though, because what's actually called
# in this test depends on the order of entries in the PATH and on whether the
# selected executable is from a wheel build or a regular OGS build: in a wheel
# build the ogs executable is a Python script - see ogs_with_args() in Applications/Python/ogs/_internal/provide_ogs_cli_tools_via_wheel.py - in a
# regular OGS build it is a standalone binary.
@pytest.mark.parametrize("ogs_use_path", ["0", "1"])
@pytest.mark.parametrize("exe", ["ogs", "generateStructuredMesh"])
def test_ogs_use_path(ogs_use_path, exe):
    print(f"\n==== test case OGS_USE_PATH={ogs_use_path} {exe}")
    env = os.environ.copy()
    env["OGS_USE_PATH"] = ogs_use_path

    # Unix needs special handling to properly terminate any subsubprocesses
    # spawned. API is system dependent, therefore two branches.
    if os.name == "posix":
        p = subprocess.Popen([exe, "--version"], env=env, start_new_session=True)
        try:
            p.wait(timeout=1)
        except subprocess.TimeoutExpired:
            os.killpg(os.getpgid(p.pid), signal.SIGTERM)
            raise

        return_code = p.poll()
    else:
        return_code = subprocess.call([exe, "--version"], timeout=1)

    assert return_code == 0
