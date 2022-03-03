Chemistry Library
============

This library is dedicated to providing an interface for interacting with the geochemical solver Phreeqc, and shall be used together with the ComponentTransport Module for reactive transport modeling in saturated porous media. At present, this interface supports chemical calculations with a variety of reaction types including kinetically/equilibrium-controlled dissolution-precipitation reactions, surface complexation reactions, and ion-exchange reactions.

If you would like to download OpenGeoSys executables already built for different platforms, please visit the [OpenGeoSys Release Page](https://www.opengeosys.org/releases/).

Before building your own OpenGeoSys model and start the simulation, please first have a look at our [User Guide](https://www.opengeosys.org/docs/userguide/basics/introduction/). There are tutorials available regarding how to download and build the source code, as well as how to run a model simulation. For advanced users who would like to make changes in the source code, the [Developer Guide](https://www.opengeosys.org/docs/devguide/getting-started/introduction/) is also very helpful.

There are a number of [benchmarks](https://github.com/ufz/ogs/tree/master/Tests/Data/Parabolic/ComponentTransport/ReactiveTransport) related to reactive transport process. Online documentation is also available for most of them. Most relevant benchmarks include
- [1D dolomitization process in a calcite-containing porous column](https://www.opengeosys.org/docs/benchmarks/reactive-transport/calcite/)
- [1D U(VI) migration with simultaneous sorption/desoprtion process](https://www.opengeosys.org/docs/benchmarks/reactive-transport/radionuclide/radionuclide/)

If the reader would like to know details of the input parameters, please also visit our [Doxygen page](https://doxygen.opengeosys.org/).

Further notes for each header file and implementation file:

- CreateChemicalSolverInterface.h and CreateChemicalSolverInterface.cpp --- For instantiating a chemical solver interface

- ChemicalSolverInterface.h --- For providing a polymorphic interface to different chemical solvers

- ChemicalSolverType.h --- For defining the chemical solvers that OpenGeoSys-6 supports

- PhreeqcIO.h and PhreeqcIO.cpp --- For preparing Phreeqc input file to and reading Phreeqc output file

- CMakeLists.txt --- For building this library

- PhreeqcIOData --- Holding Phreeqc metadata

The following header files and implementation files are used to perform reactive transport simulations with using the direct memory access approach:

- PhreeqcKernelData

- PhreeqcKernel.cpp

- PhreeqcKernel.h
