# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

import json
import subprocess
import sys

# ruff: noqa: E402
import pytest

ogstools = pytest.importorskip("ogstools")

import ogs
import ogstools as ot


def run(data, tmp_path):
    jsonfile = tmp_path / "data.json"
    with jsonfile.open("w") as f:
        json.dump(data, f, indent=4)
    bulk_mesh = tmp_path / "rect_3x3.vtu"
    anchor_mesh = tmp_path / "test.vtu"
    return subprocess.run(
        [
            "CreateAnchors",
            "-i",
            str(bulk_mesh),
            "-o",
            str(anchor_mesh),
            "-f",
            str(jsonfile),
            "--tolerance",
            "1e-6",
            "--max-iter",
            "10",
        ],
        capture_output=True,
        check=False,
    )


@pytest.mark.skipif(
    ogs.OGS_USE_PETSC == "ON",
    reason="used VTK version is likely too old",
)
@pytest.mark.skipif(
    ogs.OGS_COVERAGE == "ON",
    reason="used VTK version is likely too old",
)
@pytest.mark.skipif(
    sys.platform == "win32",
    reason="remove skip if ogstools version is updated to 0.7.2",
)
def test_error_msg(tmp_path):
    bulk_mesh_gmsh = tmp_path / "rect.msh"
    ot.meshlib.rect(
        lengths=1.0,
        n_edge_cells=3,
        n_layers=1,
        structured_grid=True,
        order=1,
        mixed_elements=False,
        jiggle=0.0,
        out_name=bulk_mesh_gmsh,
        msh_version=None,
        layer_ids=None,
    )

    ms = ot.meshlib.meshes_from_gmsh(bulk_mesh_gmsh, dim=2, reindex=True, log=True)
    ms["domain"].save(tmp_path / "rect_3x3.vtu")

    data = {
        "anchor_start_points": [[0.1, 0.2, 0.0]],
        "anchor_end_points": [[0.8, 0.7, 0.0]],
        "maximum_anchor_stress": [500e26],
        "initial_anchor_stress": [5e7],
        "residual_anchor_stress": [250e6],
        "anchor_cross_sectional_area": [0.12566370614359174],
        "anchor_stiffness": [100e9],
        "free_fraction": [0.0],
    }

    out = run(data, tmp_path)
    assert out.returncode == 0

    data2 = data.copy()
    data2.pop("free_fraction")
    out = run(data2, tmp_path)
    assert out.returncode == 0

    data2.pop("residual_anchor_stress")
    out = run(data2, tmp_path)
    assert out.returncode == 0

    data2.pop("initial_anchor_stress")
    out = run(data2, tmp_path)
    assert out.returncode == 0

    data2.pop("maximum_anchor_stress")
    out = run(data2, tmp_path)
    assert out.returncode == 0

    data2.pop("anchor_stiffness")
    out = run(data2, tmp_path)
    assert out.returncode != 0

    data2 = data.copy()
    data2.pop("anchor_cross_sectional_area")
    out = run(data2, tmp_path)
    assert out.returncode != 0

    data2 = data.copy()
    data2.pop("anchor_end_points")
    out = run(data2, tmp_path)
    assert out.returncode != 0

    data2 = data.copy()
    data2.pop("anchor_start_points")
    out = run(data2, tmp_path)
    assert out.returncode != 0

    data["anchor_start_points"] = [-3.0, 1, 0.0]
    out = run(data, tmp_path)
    assert out.returncode != 0

    data["anchor_start_points"] = [0.1, 0.2, 0.1]
    data["anchor_end_points"] = [0.3, 0.5, 0.1]
    out = run(data, tmp_path)
    assert out.returncode != 0
