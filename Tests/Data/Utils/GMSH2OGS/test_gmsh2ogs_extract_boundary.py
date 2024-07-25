import os
from pathlib import Path

import pytest
from ogs import cli


@pytest.mark.parametrize(
    ("mesh", "gmsh2_physical_id"), [("A2-gmsh.msh", True), ("A2-gmsh4.msh", False)]
)
def test_gmsh2ogs_extract_boundary(tmp_path, mesh, gmsh2_physical_id):
    os.chdir(Path(__file__).resolve().parent)
    mesh_file = Path(f"{tmp_path}/A2.vtu")

    cli.GMSH2OGS(
        gmsh2_physical_id=gmsh2_physical_id,
        i=mesh,
        o=mesh_file,
        e=True,
        boundaries=True,
        validation=True,
    )
    assert mesh_file.exists()

    # geometry bulk
    assert cli.vtkdiff(mesh_file, "A2.vtu", mesh_check=None) == 0
    assert (
        cli.vtkdiff(
            mesh_file,
            "A2.vtu",
            a="MaterialIDs",
            b="MaterialIDs",
            abs=0,
            rel=0,
        )
        == 0
    )

    for x in range(8):
        # geometry boundary
        assert (
            cli.vtkdiff(f"{tmp_path}/A2_{x}.vtu", f"A2_{x}.vtu", mesh_check=None) == 0
        )
        assert (
            cli.vtkdiff(
                f"{tmp_path}/A2_{x}.vtu",
                f"A2_{x}.vtu",
                a="bulk_node_ids",
                b="bulk_node_ids",
                abs=0,
                rel=0,
            )
            == 0
        )
        if x in [1, 2]:  # inner boundarie
            assert (
                cli.vtkdiff(
                    f"{tmp_path}/A2_{x}.vtu",
                    f"A2_{x}.vtu",
                    a="number_bulk_elements",
                    b="number_bulk_elements",
                    abs=0,
                    rel=0,
                )
                == 0
            )

        if x in [0, 3, 4, 5, 6, 7]:  # outer boundarie
            assert (
                cli.vtkdiff(
                    f"{tmp_path}/A2_{x}.vtu",
                    f"A2_{x}.vtu",
                    a="bulk_element_ids",
                    b="bulk_element_ids",
                    abs=0,
                    rel=0,
                )
                == 0
            )
