from skbuild import setup
from setuptools import find_packages

setup(
    name="OpenGeoSys",
    version="6.4.2",
    description="OpenGeoSys",
    author="OpenGeoSys Community",
    license="BSD-3-Clause",
    packages=find_packages(where="Applications/Python"),
    package_dir={"": "Applications/Python"},
    cmake_install_dir="Applications/Python/OpenGeoSys",
    extras_require={"test": ["pytest"]},
    cmake_args=["--preset wheel"],
    python_requires=">=3.6",
)
