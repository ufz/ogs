from skbuild import setup
from setuptools import find_packages

import os
import re
import subprocess
import sys


def get_version():
    git_version = ""
    if "OGS_VERSION" in os.environ:
        git_version = os.environ["OGS_VERSION"]
    else:
        git_version = subprocess.run(
            ["git", "describe", "--tags"],
            capture_output=True,
            text=True,
            shell=True,
        ).stdout.strip()

    if re.match("\d+\.\d+\.\d+-\d+-g\w+", git_version):
        # Make it PEP 440 compliant
        # e.g. 6.4.2-1140-g85bbc8b4e1 -> 6.4.2.dev1140
        m = re.match(".+?(?=-g[\w]*$)", git_version)  # strip out commit hash
        if m:
            return m.group(0).replace("-", ".dev")  # insert dev
        else:
            print("WARNING: Could not get ogs version!")
            exit(1)
    else:
        return git_version


sys.path.append(os.path.join("Applications", "Python"))
from OpenGeoSys import binaries_list

console_scripts = []
for b in binaries_list:
    console_scripts.append(f"{b}=OpenGeoSys:{b}")

import platform
cmake_preset = "wheel"
if platform.system() == "Windows":
    cmake_preset += "-win"

setup(
    name="OpenGeoSys",
    version=get_version(),
    description="OpenGeoSys",
    author="OpenGeoSys Community",
    license="BSD-3-Clause",
    packages=find_packages(where="Applications/Python"),
    package_dir={"": "Applications/Python"},
    cmake_install_dir="Applications/Python/OpenGeoSys",
    extras_require={"test": ["pytest"]},
    cmake_args=[f"--preset {cmake_preset}", "-B ."],
    python_requires=">=3.6",
    entry_points={"console_scripts": console_scripts},
)
