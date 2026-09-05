# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

from pathlib import Path

from binder import (
    BINDER_SPARSE_PATHS,
    sparse_path_query,
    sparse_paths_for_notebook,
)


def test_notebook_directory_is_the_default_sparse_path():
    notebook = Path("Tests/Data/Parabolic/Richards/richards-flow.py")

    assert sparse_paths_for_notebook(notebook) == (notebook.parent,)
    assert sparse_path_query(sparse_paths_for_notebook(notebook)) == (
        "%26sparsePath=Tests/Data/Parabolic/Richards"
    )


def test_external_notebook_dependencies_are_repeated_sparse_paths():
    notebook = Path(
        "Tests/Data/Parabolic/T/3D_line_source_term_tests/"
        "3D_line_source_term_in_cylinder/heatconduction-line_source_term.py"
    )
    paths = sparse_paths_for_notebook(notebook)

    assert paths == (
        notebook.parent,
        Path(
            "Tests/Data/Parabolic/T/3D_line_source_term_tests/"
            "3D_line_source_term_in_cylinder_axisymmetric"
        ),
    )
    assert sparse_path_query(paths) == (
        "%26sparsePath=Tests/Data/Parabolic/T/3D_line_source_term_tests/"
        "3D_line_source_term_in_cylinder"
        "%26sparsePath=Tests/Data/Parabolic/T/3D_line_source_term_tests/"
        "3D_line_source_term_in_cylinder_axisymmetric"
    )


def test_explicit_additional_paths_replace_fallback_paths():
    notebook = Path("Tests/Data/example/example.py")
    paths = sparse_paths_for_notebook(
        notebook, ("Tests/Data/shared", "Tests/Data/shared")
    )

    assert paths == (
        notebook.parent,
        Path("Tests/Data/shared"),
    )


def test_all_audited_external_notebooks_have_sparse_path_overrides():
    assert set(BINDER_SPARSE_PATHS) == {
        Path("Tests/Data/HMPhaseField/GreatCell/GreatCellHM_VPF.py"),
        Path("Tests/Data/HMPhaseField/GreatCell/GreatCellHM_VPF_propagating.py"),
        Path("Tests/Data/LIE/HydroMechanics/GreatCelljupyterNotebook/GreatCellHM.py"),
        Path(
            "Tests/Data/Parabolic/ComponentTransport/"
            "DiffusionSorptionDecay/DiffusionSorptionDecay.py"
        ),
        Path("Tests/Data/Parabolic/ComponentTransport/elder_jupyter/elder_jupyter.py"),
        Path(
            "Tests/Data/Parabolic/T/3D_line_source_term_tests/"
            "3D_line_source_term_in_cylinder/heatconduction-line_source_term.py"
        ),
    }
