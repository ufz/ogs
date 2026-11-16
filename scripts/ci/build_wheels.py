"""Build wheels using the version computed by the CI metadata job."""

import os
import re
import subprocess


def python_version(ogs_version):
    match = re.fullmatch(
        r"v?(\d+)\.(\d+)\.(\d+)(?:-?rc(\d+))?(?:-(\d+)-g[0-9a-f]+)?(-dirty)?",
        ogs_version,
    )
    if not match:
        message = f"Unsupported OGS_VERSION: {ogs_version}"
        raise ValueError(message)

    major, minor, patch, rc, distance, dirty = match.groups()
    development = bool(int(distance or 0) or dirty)
    version = f"{major}.{minor}.{patch}"
    if rc is not None:
        version += f"rc{int(rc) + development}"
    elif development:
        version = f"{major}.{minor}.{int(patch) + 1}"
    if development:
        version += f".dev{int(distance or 0)}"
    return version


if __name__ == "__main__":
    environment = os.environ.copy()
    environment["SETUPTOOLS_SCM_PRETEND_VERSION_FOR_OGS"] = python_version(
        environment["OGS_VERSION"]
    )
    print(
        "Python package version:",
        environment["SETUPTOOLS_SCM_PRETEND_VERSION_FOR_OGS"],
        flush=True,
    )
    subprocess.run(["uvx", "cibuildwheel"], env=environment, check=True)
