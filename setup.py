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
    cmake_args=[
        "-DOGS_BUILD_PROCESSES=LiquidFlow",
        "-DOGS_BUILD_UTILS=OFF",
        "-DHDF5_USE_STATIC_LIBRARIES=ON",
        "-DOGS_BUILD_HDF5=ON",
        "-DOGS_USE_PYTHON=OFF",  # not possible because manylinux does not provide libpythonX.Y.so
        "-DOGS_BUILD_PYTHON_MODULE=ON",
        "-DOGS_BUILD_TESTING=OFF",
        "-DOGS_INSTALL_DEPENDENCIES=OFF",  # otherwise auditwheel fails
    ],
    python_requires=">=3.6",
)
