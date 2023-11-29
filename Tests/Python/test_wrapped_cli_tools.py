import os
import tempfile

import ogs


def test_ogs_version():
    assert ogs.cli.ogs("--version") == 0
    assert ogs.cli.ogs(version=None) == 0


def test_ogs_help():
    assert ogs.cli.ogs("--help") == 0
    assert ogs.cli.ogs(help=None) == 0


def test_generate_structured_mesh():
    with tempfile.TemporaryDirectory() as tmpdirname:
        outfile = os.path.join(tmpdirname, "test.vtu")
        assert not os.path.exists(outfile)

        assert ogs.cli.generateStructuredMesh(e="line", lx=1, nx=10, o=outfile) == 0

        assert os.path.exists(outfile)


def test_parameter_with_underscore():
    tempdir = tempfile.mkdtemp()
    gml_filename = os.path.join(tempdir, "geom.gml")
    gmsh_filename = os.path.join(tempdir, "geom.gmsh")

    assert (
        ogs.cli.generateGeometry(
            geometry_name="SceneRectangle",
            x0="0",
            x1="10",
            y0="0",
            y1="10",
            z0="0",
            z1="0",
            nx="0",
            nx1="0",
            ny="0",
            ny1="0",
            nz="0",
            nz1="0",
            polyline_name="bounding_rectangle",
            o=gml_filename,
        )
        == 0
    )
    assert os.path.exists(gml_filename)

    assert (
        ogs.cli.geometryToGmshGeo(
            i=gml_filename, o=gmsh_filename, mesh_density_scaling_at_points=20
        )
        == 0
    )
    assert os.path.exists(gmsh_filename)
