import json
import subprocess

# ruff: noqa: E402
import pytest

ogstools = pytest.importorskip("ogstools")

import ogs
import ogstools as ot


@pytest.mark.skipif(
    ogs.OGS_USE_PETSC == "ON",
    reason="used VTK version is likely too old",
)
@pytest.mark.skipif(
    ogs.OGS_COVERAGE == "ON",
    reason="used VTK version is likely too old",
)
@pytest.mark.parametrize(
    "order_intersections", [(1, [6, 6, 6, 6, 10]), (2, [6, 6, 6, 6, 10])]
)
def test_number_intersection(tmp_path, order_intersections):
    order, intersections = order_intersections
    bulk_mesh_gmsh = tmp_path / "rect.msh"
    ot.meshlib.rect(
        lengths=1.0,
        n_edge_cells=3,
        n_layers=1,
        structured_grid=True,
        order=order,
        mixed_elements=False,
        jiggle=0.0,
        out_name=bulk_mesh_gmsh,
        msh_version=None,
        layer_ids=None,
    )

    ms = ot.meshlib.meshes_from_gmsh(bulk_mesh_gmsh, dim=2, reindex=True, log=True)
    bulk_mesh = tmp_path / f"rect_3x3_{order}.vtu"
    ms["domain"].save(bulk_mesh)

    startp = [
        [1.0e-15, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.12, 0.3, 0.0],
    ]
    endp = [
        [1.0e-15, 1.0, 0.0],
        [1.0, 1.0, 0.0],
        [1.0 + 1e-15, 1.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.8, 0.7, 0.0],
    ]

    for i, entry in enumerate(startp):
        data = {
            "anchor_start_points": [entry],
            "anchor_end_points": [endp[i]],
            "maximum_anchor_stress": [500e26],
            "initial_anchor_stress": [5e7],
            "residual_anchor_stress": [250e6],
            "anchor_cross_sectional_area": [0.12566370614359174],
            "anchor_stiffness": [100e9],
            "free_fraction": [0.0],
        }
        json_file = tmp_path / "data.json"
        with json_file.open("w") as f:
            json.dump(data, f, indent=4)
        anchor_mesh = tmp_path / f"test_{order}.vtu"
        subprocess.run(
            [
                "CreateAnchors",
                "-i",
                str(bulk_mesh),
                "-o",
                str(anchor_mesh),
                "-f",
                str(json_file),
                "--tolerance",
                "1e-12",
                "--max-iter",
                "10",
            ],
            capture_output=True,
            check=True,
        )
        anchor_mesh = ot.Mesh(anchor_mesh)
        assert intersections[i] == anchor_mesh.n_points
        assert intersections[i] == 2 * anchor_mesh.n_cells
