# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

from pathlib import Path
from urllib.parse import quote

# Additional repository paths needed by notebooks when Binder performs a sparse
# checkout. The notebook directory itself is added automatically below.
BINDER_SPARSE_PATHS: dict[Path, tuple[Path, ...]] = {
    Path("Tests/Data/HMPhaseField/GreatCell/GreatCellHM_VPF.py"): (
        Path("Tests/Data/LIE/Mechanics/GreatCelljupyterNotebook"),
    ),
    Path("Tests/Data/HMPhaseField/GreatCell/GreatCellHM_VPF_propagating.py"): (
        Path("Tests/Data/LIE/Mechanics/GreatCelljupyterNotebook"),
    ),
    Path("Tests/Data/LIE/HydroMechanics/GreatCelljupyterNotebook/GreatCellHM.py"): (
        Path("Tests/Data/LIE/Mechanics/GreatCelljupyterNotebook"),
    ),
    Path(
        "Tests/Data/Parabolic/ComponentTransport/"
        "DiffusionSorptionDecay/DiffusionSorptionDecay.py"
    ): (
        Path("Tests/Data/Parabolic/ComponentTransport/AdvectionDiffusionSorptionDecay"),
    ),
    Path("Tests/Data/Parabolic/ComponentTransport/elder_jupyter/elder_jupyter.py"): (
        Path("Tests/Data/Parabolic/ComponentTransport/elder"),
    ),
    Path(
        "Tests/Data/Parabolic/T/3D_line_source_term_tests/"
        "3D_line_source_term_in_cylinder/heatconduction-line_source_term.py"
    ): (
        Path(
            "Tests/Data/Parabolic/T/3D_line_source_term_tests/"
            "3D_line_source_term_in_cylinder_axisymmetric"
        ),
    ),
}


def sparse_paths_for_notebook(
    notebook: Path, additional_paths: tuple[str | Path, ...] = ()
) -> tuple[Path, ...]:
    """Return the notebook directory and any additional Binder paths."""
    if additional_paths:
        extra_paths = tuple(Path(path) for path in additional_paths)
    else:
        # Keep direct testrunner invocations working until all callers provide
        # the dependencies explicitly through NotebookTest().
        extra_paths = BINDER_SPARSE_PATHS.get(notebook, ())

    paths = (notebook.parent, *extra_paths)
    return tuple(dict.fromkeys(paths))


def sparse_path_query(paths: tuple[Path, ...]) -> str:
    """Build repeated, encoded sparsePath parameters for the inner URL."""
    return "".join(
        f"%26sparsePath={quote(path.as_posix(), safe='/')}" for path in paths
    )
