#!/usr/bin/env python

import os

# need to increase OGS_CTEST_MAX_RUNTIME to enable these:
# ctests = ["SurfaceComplexation", "EquilibriumPhase", "KineticReactant"]
ctests = ["SteadyState", "ComponentTransport", "ThermoHydroMechanics"]
report_path = "./coverage_reports"

os.makedirs(report_path, exist_ok=True)

for t in ctests:
    os.system("cmake --build . -t clean_coverage")
    os.system(f"ctest -R {t}")
    report = f"{report_path}/{t}_coverage.json"
    os.system(
        f"@FASTCOV_PATH@ --branch-coverage --include @PROJECT_SOURCE_DIR@ --gcov @GCOV_PATH@ --search-directory @PROJECT_BINARY_DIR@ --process-gcno --output {report} --exclude Applications/CLI/ ProcessLib/ Tests/ --exclude Applications/CLI/ ProcessLib/ Tests/"
    )
    os.system(f"perl -i -pe s!@PROJECT_SOURCE_DIR@/!!g {report}")
    if "CI" in os.environ:
        print(
            f"Coverage visualization: https://ogs.ogs.xyz/code-coverage-vis/?report={os.environ.get('CI_JOB_URL')}/artifacts/raw/build/coverage/{report}"
        )
