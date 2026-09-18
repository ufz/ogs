# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

import os
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.skipif(sys.platform != "win32", reason="Windows DLL loading")
@pytest.mark.skipif("OGS_USE_PATH" in os.environ, reason="Works in wheel only.")
@pytest.mark.parametrize("module", ["callbacks", "mpl", "OGSMesh", "OGSSimulator"])
def test_extension_import_without_build_environment(module, tmp_path):
    # A fresh interpreter must find the bundled DLLs without the test suite's
    # oneAPI directories, previously loaded extensions, or the build PATH.
    env = {
        key: value
        for key, value in os.environ.items()
        if key.upper() not in {"PATH", "PYTHONPATH", "PYTHONHOME", "MKLROOT"}
    }
    env["PATH"] = str(Path(os.environ["SYSTEMROOT"]) / "System32")
    subprocess.run(
        [sys.executable, "-I", "-c", f"import ogs.{module}"],
        cwd=tmp_path,
        env=env,
        check=True,
    )
