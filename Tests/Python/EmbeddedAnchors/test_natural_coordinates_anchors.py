import json
import pathlib
import random
import subprocess

# ruff: noqa: E402
import pytest

ogstools = pytest.importorskip("ogstools")

import numpy as np
import ogstools as ot


@pytest.mark.parametrize("order", [1, 2])
def test_natural_coordinates(tmp_path, order):
    ot.meshlib.rect(
        lengths=1.0,
        n_edge_cells=2,
        n_layers=1,
        structured_grid=True,
        order=order,
        mixed_elements=False,
        jiggle=0.0,
        out_name=tmp_path / pathlib.Path("rect.msh"),
        msh_version=None,
        layer_ids=None,
    )

    ms = ot.meshlib.meshes_from_gmsh(
        tmp_path / "rect.msh", dim=2, reindex=True, log=True
    )
    ms["domain"].save(tmp_path / "rect_2x2.vtu")

    for _ in range(10):
        x1 = random.random()
        y1 = random.random()
        x2 = random.random()
        y2 = random.random()

        data = {
            "anchor_start_points": [[x1, y1, 0]],
            "anchor_end_points": [[x2, y2, 0]],
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

        bulk_mesh = tmp_path / "rect_2x2.vtu"
        anchor_mesh = tmp_path / "test.vtu"
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
                "100",
            ],
            capture_output=True,
            check=True,
        )

        anchor_mesh = ot.Mesh(anchor_mesh)
        cell_data = [
            "initial_anchor_stress",
            "maximum_anchor_stress",
            "residual_anchor_stress",
            "anchor_cross_sectional_area",
            "anchor_stiffness",
        ]
        point_data = [
            "natural_coordinates",
            "bulk_element_ids",
            "point_cloud_node_ids",
        ]
        for entry in cell_data:
            assert entry in anchor_mesh.cell_data

        n_cells = len(anchor_mesh.cell_data["anchor_stiffness"])

        for entry in point_data:
            assert entry in anchor_mesh.point_data
        n_points = len(anchor_mesh.points)
        assert n_points == 2 * n_cells
        subtract = {
            0: [0.25, 0.25],
            1: [0.25, 0.75],
            2: [0.75, 0.25],
            3: [0.75, 0.75],
        }
        delta = 1e-11
        for i, point in enumerate(anchor_mesh.points):
            px = (
                -(point[0] - subtract[anchor_mesh.point_data["bulk_element_ids"][i]][0])
                * 4
            )
            py = (
                -(point[1] - subtract[anchor_mesh.point_data["bulk_element_ids"][i]][1])
                * 4
            )
            assert (
                np.abs(anchor_mesh.point_data["natural_coordinates"][i][0] - px) < delta
            )
            assert (
                np.abs(anchor_mesh.point_data["natural_coordinates"][i][1] - py) < delta
            )
