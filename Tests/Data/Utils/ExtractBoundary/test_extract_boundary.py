import os
from pathlib import Path

import pytest
from ogs import cli


@pytest.mark.parametrize(("element_type"), [("tri"), ("quad")])
def test_extract_boundary(tmp_path, element_type):
    os.chdir(Path(__file__).resolve().parent)
    mesh_basename = f"square_10_1x1_{element_type}"
    mesh_file = Path(f"{tmp_path}/{mesh_basename}.vtu")
    cli.generateStructuredMesh(
        e=element_type,
        lx=1,
        ly=1,
        nx=10,
        ny=10,
        o=mesh_file,
    )
    assert mesh_file.exists()

    boundary_file = Path(f"{tmp_path}/{mesh_basename}_boundary.vtu")
    cli.ExtractBoundary(i=mesh_file, o=boundary_file)
    assert boundary_file.exists()

    assert cli.vtkdiff(boundary_file, boundary_file.name, mesh_check=True) == 0
    for field in ["bulk_node_ids", "bulk_element_ids", "bulk_face_ids"]:
        assert (
            cli.vtkdiff(
                boundary_file,
                boundary_file.name,
                a=field,
                b=field,
                abs=0,
                rel=0,
            )
            == 0
        )
