import pytest

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

        assert 0 == ogs.cli.generateStructuredMesh(e="line", lx=1, nx=10, o=outfile)

        assert os.path.exists(outfile)
