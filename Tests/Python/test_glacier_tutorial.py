# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

import os
import platform
import runpy
import shutil
import subprocess

import ogs
import ogs._internal.provide_ogs_cli_tools_via_wheel as ogs_cli_wheel
import pytest

from . import push_argv

# Needs ogstools for msh2vtu
pytest.importorskip("ogstools")


def _run(program, args):
    func = getattr(ogs_cli_wheel, program)
    with push_argv([f"{program}.py", *args]), pytest.raises(SystemExit) as excinfo:
        func()
    assert excinfo.value.code == 0


@pytest.mark.skipif(ogs.OGS_USE_PETSC == "ON", reason="Not supported")
@pytest.mark.skipif(
    platform.system() == "Windows", reason="TODO: NodeReordering(.exe) not found"
)
def test_glacier_tutorial(tmp_path, ogs_src_dir, monkeypatch):
    monkeypatch.setenv("CI", "1")

    tutorial_dir = (
        ogs_src_dir / "web" / "content" / "docs" / "tutorials" / "advancing-glacier"
    )

    for filename in ("mesh_basin.py", "OGSinput_basin.prj", "timeBCs_glacier.py"):
        shutil.copy(tutorial_dir / filename, tmp_path / filename)

    os.chdir(tmp_path)
    runpy.run_path("mesh_basin.py", run_name="__main__")
    subprocess.run(
        ["msh2vtu", "mesh_basin.msh", "--reindex", "-o", tmp_path], check=True
    )
    _run("ogs", ["OGSinput_basin.prj"])

    assert (tmp_path / "mesh_basin_domain.vtu").exists()
    assert (tmp_path / "OGSoutput_basin.pvd").exists()
