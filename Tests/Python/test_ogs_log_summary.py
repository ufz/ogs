import filecmp
from pathlib import Path

import pandas as pd
import pandas.testing as pd_test
import pytest
from ogs.dev.ogs_log_summary import aggregate_log_files, run

res_path = Path(__file__).parent / "res/test_ogs_log_summary"


@pytest.mark.parametrize(
    ("log_file", "prj_file_expected"),
    [
        (
            "arch-ci-job-ogs-ThermoRichardsMechanics_MFront_BentoniteBehaviourGeneralMod230824_02-ramp_bentonite_ramped_Neumann_BC.txt",
            "ThermoRichardsMechanics/MFront/BentoniteBehaviourGeneralMod230824/02-ramp/bentonite_ramped_Neumann_BC.prj",
        ),
        (
            "petsc-ci-job-ogs-ThermoRichardsMechanics_MFront_BentoniteBehaviourGeneralMod230824_02-ramp_bentonite_ramped_Neumann_BC.txt",
            "ThermoRichardsMechanics/MFront/BentoniteBehaviourGeneralMod230824/02-ramp/bentonite_ramped_Neumann_BC.prj",
        ),
    ],
)
def test_aggregation(log_file, prj_file_expected):
    log_path = res_path / log_file

    map_prj_file_to_agg_vtkdiff_stats = aggregate_log_files([log_path])

    assert len(map_prj_file_to_agg_vtkdiff_stats) == 1

    assert prj_file_expected in map_prj_file_to_agg_vtkdiff_stats

    df_agg_vtkdiff_stats = map_prj_file_to_agg_vtkdiff_stats[prj_file_expected]

    with Path(res_path / f"{log_file}.df_ref.csv").open() as fh:
        df_agg_vtkdiff_stats_expected = pd.read_csv(fh, sep="\t", index_col=0)

    pd_test.assert_frame_equal(
        df_agg_vtkdiff_stats, df_agg_vtkdiff_stats_expected, atol=1e-15, rtol=0
    )

    files = df_agg_vtkdiff_stats["file"].to_numpy()
    assert (
        files == files[0]
    ).all()  # the "maximum" file name must be the same for all fields


@pytest.mark.parametrize(
    ("log_file", "snippet_file_expected"),
    [
        (
            "arch-ci-job-ogs-ThermoRichardsMechanics_MFront_BentoniteBehaviourGeneralMod230824_02-ramp_bentonite_ramped_Neumann_BC.txt",
            "ThermoRichardsMechanics/MFront/BentoniteBehaviourGeneralMod230824/02-ramp/bentonite_ramped_Neumann_BC.prj",
        ),
        (
            "petsc-ci-job-ogs-ThermoRichardsMechanics_MFront_BentoniteBehaviourGeneralMod230824_02-ramp_bentonite_ramped_Neumann_BC.txt",
            "ThermoRichardsMechanics/MFront/BentoniteBehaviourGeneralMod230824/02-ramp/bentonite_ramped_Neumann_BC.prj",
        ),
    ],
)
def test_xml_snippet(tmp_path, log_file, snippet_file_expected):
    # Note: this test uses pytest's tmp_path fixture, cf. https://docs.pytest.org/en/7.1.x/how-to/tmp_path.html
    log_path = res_path / log_file

    run([log_path], tmp_path, False)

    snippet_path = tmp_path / snippet_file_expected
    assert snippet_path.exists()

    snippet_path_expected = res_path / f"{log_file}.snippet_ref.xml"
    assert filecmp.cmp(snippet_path, snippet_path_expected, shallow=False)
