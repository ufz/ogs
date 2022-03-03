Component Tranpsort Module
============

This Component Tranpsort Module is dedicated to the modeling of Darcy-scale subsurface flow and (reactive) solute transport in saturated porous media.

If you would like to download OpenGeoSys executables already built for different platforms, please visit the [OpenGeoSys Release Page](https://www.opengeosys.org/releases/).

Before building your own OpenGeoSys model and start the simulation, please first have a look at our [User Guide](https://www.opengeosys.org/docs/userguide/basics/introduction/). There are tutorials available regarding how to download and build the source code, as well as how to run a model simulation. For advanced users who would like to make changes in the source code, the [Developer Guide](https://www.opengeosys.org/docs/devguide/getting-started/introduction/) is also very helpful.

In order to ensure the reliability of the model feature, the code implementation is tested continuously through a variety of benchmarks. The input files of the benchmarks related to the Component Transport Module can be found [here](https://github.com/ufz/ogs/tree/master/Tests/Data/Parabolic/ComponentTransport). Additionally, online documentation is available for some featured benchmarks to demonstrate the capabilities of this module. For example, interested readers may refer to the [1D/2D tracer test](https://www.opengeosys.org/docs/benchmarks/hydro-component/contracer/contracer/), and [1D dolomitization process in a calcite-containing porous column](https://www.opengeosys.org/docs/benchmarks/reactive-transport/calcite/) case. If the reader would like to know details of the input parameters, please also visit our [Doxygen page](https://doxygen.opengeosys.org/).

Further notes for each header file and implementation file:

- CreateComponentTransportProcess.h and CreateComponentTransportProcess.cpp --- For instantiating solute transport process

- ComponentTransportProcessData.h --- For holding the process data that is unique to the solute transport process

- ComponentTransportProcess.h and ComponentTransportProcess.cpp --- For any process-related manipulations at the global level, e.g., process initialization, initial condition assignment, call to assemble local matrices for solving linearized PDEs, and computation of secondary variables.

- ComponentTransportFEM.h --- For any process-related manipulations at the element level

- CMakeLists.txt --- For building this process module

- Tests.cmake --- For automatically running tests in the continuous integration process

The following header files and implementation files are used to perform reactive transport simulations with using look-up table approach:
- CreateLookupTable.cpp
- CreateLookupTable.h
- LookupTable.cpp
- LookupTable.h
