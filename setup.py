from skbuild import setup
from setuptools import find_packages

import os
import platform
import sys

sys.path.append(os.path.join("Applications", "Python"))
from ogs._internal.provide_ogs_cli_tools_via_wheel import binaries_list

console_scripts = []
for b in binaries_list:
    console_scripts.append(f"{b}=ogs._internal.provide_ogs_cli_tools_via_wheel:{b}")

cmake_preset = "wheel"
if platform.system() == "Windows":
    cmake_preset += "-win"

from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

# setuptools_scm config local_scheme via env var SETUPTOOLS_SCM_LOCAL_SCHEME:
scm_local_scheme = "node-and-date"
if "SETUPTOOLS_SCM_LOCAL_SCHEME" in os.environ:
    local_scheme_values = [
        "node-and-date",
        "node-and-timestamp",
        "dirty-tag",
        "no-local-version",
    ]
    if os.environ["SETUPTOOLS_SCM_LOCAL_SCHEME"] in local_scheme_values:
        scm_local_scheme = os.environ["SETUPTOOLS_SCM_LOCAL_SCHEME"]

setup(
    name="ogs",
    description="OpenGeoSys Python Module",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="OpenGeoSys Community",
    license="BSD-3-Clause",
    packages=find_packages(where="Applications/Python"),
    package_dir={"": "Applications/Python"},
    cmake_install_dir="Applications/Python/ogs",
    extras_require={"test": ["pytest", "numpy"]},
    cmake_args=[f"--preset {cmake_preset}", "-B ."],
    python_requires=">=3.7",
    entry_points={"console_scripts": console_scripts},
    use_scm_version={
        "write_to": "Applications/Python/_version.py",
        "write_to_template": '__version__ = "{version}"',
        # Switching to guess-next-dev (default) would increment the version number
        # This would be in line with PEP 440, switch OGS versioning too?
        "version_scheme": "no-guess-dev",
        "local_scheme": scm_local_scheme,
    },
)
