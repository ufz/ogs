# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

import os
from pathlib import Path

from ogs import cli


def test_voxelgrid_from_layers(tmp_path):
    os.chdir(Path(__file__).resolve().parent)

    mesh_file = Path(f"{tmp_path}/AREHS_test.vtu")

    cli.Layers2Grid(i="AREHS_test_layers.txt", o=mesh_file, x=500, y=300, z=100)
    assert mesh_file.exists()

    assert cli.vtkdiff(mesh_file, "AREHS_test.vtu", mesh_check=None) == 0

    fault_mesh_file = Path(f"{tmp_path}/AREHS_test_fault.vtu")
    cli.AddFaultToVoxelGrid(i=mesh_file, f="AREHS_fault.vtu", o=fault_mesh_file)
    assert fault_mesh_file.exists()

    assert cli.vtkdiff(fault_mesh_file, "AREHS_test_fault.vtu", mesh_check=None) == 0
    assert (
        cli.vtkdiff(
            fault_mesh_file,
            "AREHS_test_fault.vtu",
            a="MaterialIDs",
            b="MaterialIDs",
            abs=0,
            rel=0,
        )
        == 0
    )

    # iso
    mesh_file = Path(f"{tmp_path}/AREHS_test_iso.vtu")

    cli.Layers2Grid(i="AREHS_test_layers.txt", o=mesh_file, x=500)
    assert mesh_file.exists()

    assert cli.vtkdiff(mesh_file, "AREHS_test_iso.vtu", mesh_check=None) == 0
