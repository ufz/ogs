import os
from pathlib import Path

import pytest


def pytest_runtest_setup(item):
    perf_test_markers = list(item.iter_markers("performance_test"))
    skip_perf_tests = os.getenv(
        "OGS_PERFORMANCE_TESTS_ALLOWED_TO_FAIL", "False"
    ).lower() in (
        "true",
        "1",
        "t",
    )

    if perf_test_markers and skip_perf_tests:
        item.add_marker(
            pytest.mark.xfail(
                reason="the environment variable OGS_PERFORMANCE_TESTS_ALLOWED_TO_FAIL says that performance tests might fail",
                run=True,
                strict=False,
            )
        )


@pytest.fixture(scope="session")
def ogs_src_dir():
    return Path(__file__).parent.parent.parent
