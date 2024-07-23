import itertools
import os
import tempfile
from pathlib import Path

from ogs import cli


def test_partmesh_mixed_elements():
    os.chdir(Path(__file__).resolve().parent)

    meshes = ["linear", "quadratic"]
    partitions = [2, 4, 8, 12]

    for mesh, partition in itertools.product(meshes, partitions):
        with tempfile.TemporaryDirectory() as tmpdirname:
            cli.partmesh(ogs2metis=None, i=f"{mesh}_mesh.vtu", o=tmpdirname)
            assert Path(f"{tmpdirname}/{mesh}_mesh.mesh").exists()

            cli.partmesh(
                exe_metis=None,
                i=f"{mesh}_mesh.vtu",
                x=f"{tmpdirname}/{mesh}_mesh",
                o=tmpdirname,
                np=partition,
            )
            for filetype in [
                "cell_properties_cfg",
                "cell_properties_val",
                "msh_cfg",
                "msh_ele_g",
                "msh_ele",
                "msh_nod",
            ]:
                assert Path(
                    f"{tmpdirname}/{mesh}_mesh_partitioned_{filetype}{partition}.bin"
                ).exists()
