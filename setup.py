from skbuild import setup
from setuptools import find_packages

setup(
    name="ogs",
    version="6.4.2",
    description="OpenGeoSys",
    author="OpenGeoSys Community",
    license="BSD-3-Clause",
    packages=find_packages(where="."),
    package_dir={"": "."},
    cmake_args=["-DOGS_BUILD_PROCESSES=SteadyStateDiffusion", "-DOGS_BUILD_UTILS=OFF"],
    python_requires=">=3.6",
)
