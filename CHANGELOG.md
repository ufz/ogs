## Release notes

# 6.0.9 (In preparation)

### Features

### Utilities

### Infrastructure

- CMake option OGS_EIGEN_DYNAMIC_SHAPE_MATRICES defaults to OFF on Release
  config, ON otherwise. Can be overridden by explicitly setting the option. #1673

### Fixes

# 6.0.8

The highlight of the release is the implementation of the Lower-Interface
Elements for both the small deformation process (M) and hydro-mechanics process
(HM) allowing fractures to be incorporated in the solution domain.
For the liquid flow and two-phase flow processes several material models for the
fluids pressure, density, permeability, and viscosity were added.

### Features:

 - Implementation of hydro-mechanics (HM) with LIE. #1537-#1541
 - Implementation of small deformation (M) with LIE. #1452
 - Fracture constitutive models. #1434
 - Hydro-Mechanics process. #1508
 - First version of monolithic hydro-thermal process implementation with
   Boussinesq approximation using constant viscosity. #1534
 - Two phase flow process with pp model implementation. #1530
 - Richards flow process implementation. #1473
 - Liquid process. #1468
 - Classes for relative permeability models. #1531
 - Classes for capillary models. #1517, #1578
 - Ehlers single-surface yield function constitutive relation model. #1556
 - Support scaling, GMRES, and Pardiso in Eigen linear solvers. #1509 #1510
 - Piecewise linear Monotonic curve and a generic curve parser. #1529
 - Support searching boundary nodes in MeshLib::NodeSearch. #1459
 - Support specifying the shape function order in process variables. #1503
 - Command line option --unbuffered-std-out to deactivate buffer for standard output. #1514
 - CMake option OGS_FATAL_ABORT for debugging. #1432
 - Set the default OGS_LOG_LEVEL to debug in release builds. #1522
 - Add integration order in input files. #1464

### Utilities
New utilities:
 - createQuadraticMesh #1500
 - convertToLinearMesh #1554
 - postLIE #1555

New features:
 - extend NodeReordering to correct ordering of nonlinear nodes #1519


### Infrastructure:

 - Ctest now works on Windows too by removing time-wrappers. #1480
 - Moved to public Jenkins instance at jenkins.opengeosys.org. #1505
 - Doxygen warnings parser in Jenkins will mark a build as unstable
   if there are Doxygen warnings. #1585
 - Benchmarking on Jenkins now saves the standard output into a file for each
   test. #1528

### Fixes:
 - Fix LocalToGlobalIndexMap with mutliple variables and with multiple componets. #1433 #1440
 - Fix PropertyVector<T*> for multi-component case. #1441
 - Fix checking IDs of nonlinear nodes. #1495
 - Fix incorrect use of getNumberOfBaseNodes(). #1515
 - Fix computing sparsity pattern for mixed shape function order cases. #1548
 - Fix that iterations and residuals were not printed when Eigen linear solver fails. #1499
 - Fix all of the Doxygen warnings in the code. #1569 #1573
 - Fix all of the input file/keyword documentation and its generation.


# 6.0.7

### Features:

The main features of this release is the implementation of two new processes,
the small deformation, and the heat conduction.  Some extensions were done to
the DOF table to be able to manage multi-component/multi-variable processes.
Also, during implementation of the Robin boundary conditions, the base classes
for boundary conditions were generalized.

 - Add small deformation process with linear elastic material model. The
   implementation is based on the Kelvin mapping. #1340
 - Added B-Matrices and Kelvin mapping tools for deformation processes. #1359
 - Heat conduction process implementation. #1328
 - Finalize support for multicomponent boundary conditions adding configuration
   parser. #1343
 - Added uniform Robin boundary condition. #1336
 - Added a generic natural boundary condition class. #1337
 - Added Robin boundary condition. #1336
 - Reworked the Parameter class. It now serves as a basis for BCs and ICs.
   #1357, #1356
 - Added time-dependent Dirichlet BCs. #1380
 - Add calculation of surface flux, tests for groundwater flow. #1429
 - Implemented numerical Jacobian assembly for Newton-Raphson solver. #1352
 - Added the new parameter type "Group" which can be used for setting material
   ID dependent values. #1426
 - Added fluid property class and several fluid density and viscosity models
   based on it. #1398, #1435
 - Enabled solving of axially symmetric problems on 2D meshes for all currently
   implemented processes. #1443
 - Added time measurement for assembly and solvers. #1322
 - Added named functions, out of which expressions can be built up at run-time
   from the prj file, which can be used to output additional nodal quantities.
   #1314, #1315
 - Added component-wise norms, and flexible convergence criteria for nonlinear
   solvers. #1349, #1342
 - Restructured the time loop. #1364

### Utilities
New utilities:
 - createNeumannBc: The tool integrates the given element property and writes
   the computed data as a PropertyVector with the name `node_aggregated_gwn`
   into the mesh. The tool also outputs an OGS-5 direct source term (Neumann
   boundary condition) data file. #1346
 - scaleProperty for simple rescaling of mesh properties. #1347
 - convertGEO for geometric file conversion, e.g. gli to glm. #1360
 - swapNodeCoordinateAxes to swap node coordinate values, e.g. XY to XZ plane.
   #1361
New feature:
 - Support tetrahedra types in generateStructuredMesh. #1353

### Infrastructure:
 - Migrated all important Jenkins jobs to script-based
   [Jenkins Pipeline](https://jenkins.io/doc/pipeline/)
   functionality.  For an introduction see
   [docs.opengeosys.org - Continuous Integration](https://docs.opengeosys.org/docs/devguide/development-workflows/continuous-integration).
   #1392, #1396, #1404, #1411, #1424, #1428, #1436
 - Moved CMake logic for packaging executable dependencies (such as shared libs)
   to the install and package targets instead of running after each executable
   gets build. #1458
 - Increase minimum clang compiler version to 3.5 in course of updating travis
   build enviroment to Ubuntu LTS 14.04. #1448
 - Added a script that generates crosslinked Doxygen pages out of ctest input
   files #1348

### Fixes:
 - Fix an issue that a shape vector was defined as a column vector. Corrected to
   a row vector . #1288
 - Fix usage of `boost::optional<T const&>`, which has changed in version 1.61.
   #1385
 - Fix Grid (enlarge bounding box to fit all points). #1369
 - Fix mapping of geometries to mesh surfaces. #1327. #1368
 - Fix transmitting raster data to element properties. #1347
 - Fixed missing XSD files in packages. #1410
 - Fix a shape vector to a row vector. #1288
 - Fix FEFLOW import. #1397
 - Fix NodeReordering to check ordering of each element. #1425
 - Fix builds linking shared VTK library. #1431
 - Fix global Newton iteration counter. #1341
 - Correct few loops over mesh nodes, which should run over the mesh subsets.
   #1437
 - Fix shape function computation for 2D elements lying in the x-y-plane #1318
 - Fix AddTest, s.t. ctest now really checks results. #1325
 - Made Eigen preconditioner configurable. #1367

# 6.0.6

### Features:
 - Add external ode-solver interface with [Sundials CVODE
   library](http://computation.llnl.gov/projects/sundials-suite-nonlinear-differential-algebraic-equation-solvers/cvode). #1109
 - Add piecewise linear curves parser to the project files. The curves are
   specified by two vectors, the coordinates and values. They can be used for
   example to map temporal dependencies (time-dependent boundary conditions) or
   as approximations of coefficient dependencies (e.g. pressure-saturation
   curves). #1149
 - Extend the LocalAssemblerInterface by adding default implementations of
   pre/postTimestep and assembleJacobian functions. The current time and time
   step size are passed in the preTimestep call to the particular processes. #1214
 - Add support multi-variable/multi-component in the DOF table interface and
   extend the initial conditions to multi-components. #1224
 - Major rework of the general process interface. That also affects the
   interface of concrete processes and local assemblers.
   least squares optimization. #1145
 - Added functionality for the output of secondary variables. #1171
 - Added material properties for zeolite adsorption and Calcium oxide/hydroxide
   reactions. #1178
 - Transferred the TES process, a monolithically coupled THC model for simulating
   thermochemical energy storag devices, from OGS5. #1181
 - Introduced a general scheme for documenting OGS6 input file settings. #978
 - Added copy constructor for the class Surface, minor improvements in GeoLib. #1237
 - Added classes GeoLib::LineSegment and GeoLib::Polyline::SegmentIterator. #1139
 - GMSHInterface can handle stations as constraints. #1212
 - Added functionality to duplicate geometric data. #1229
 - Station names can be modified in Data Explorer #1273

### Infrastructure
 - Fix circular dependencies on library level. This allows for dynamic linking
   which is faster than static and can be used in debug builds, where the
   compilation time is more important than the runtime.
   - Enable shared linking of ogs libraries. #1133
   - Break FileIO on ApplicationsLib dependency. #1138
   - Remove MeshLib on FileIO dependency. #1143, #1153
   - Cleanup some of AssemblerLib dependencies. #1166
   - Split AssemblerLib and move to MathLib and NumLib #1208
   - Move InsituLib to MeshLib #1208
   - Remove MathLib depends on NumLib #1208
   - Remove dependency of FileIO on Data Explorer #1302
 - Introduced Conan package manager for automatic fetching of build dependencies, #1141
 - Inconsistent formatting of tabs and spaces was finally resolved: now all
   formatting, indentation and alignment, are done with four spaces. #1201
 - Windows 32-bit builds are disallowed because they are not supported.
   Can be forced by setting OGS_32_BIT=ON. #1230
 - Simplified FindEigen.cmake, #1209
 - git diff --check is run in its own Travis job, #1207

 - Moved some IO implementations from FileIO to BaseLib/IO, GeoLib/IO, MeshLib/IO, #1182, #1235
 - Eigen is not optional anymore #1218
 - Removed OGS_USE_EIGENLIS CMake option. Use OGS_USE_LIS instead #1251

### Fixes:
 - Fix linking of Metis in MathLib. #1147
 - Fix memory leaks in GMSHInterface. #1212
 - Fix build with Lis #1267
 - Fixing several small issues with NetCDF import #1169
 - Restructure Applications related modules
    - Move DataHolderLib and FileIO under Applications #1279
 - Remove calling std::abort within libraries. Exeptions are thrown instead. #1245
 - Fix finding Boost with BOOST_ROOT CMake argument #1287
 - Fix linking of Sundials CVODE library #1197
 - Fixed issue where geometry names would be missing after merging geometries #1295

# 6.0.5

### Features:
 - Added an ODE solver library that can solve transient and nonlinear processes
   (see http://doxygen.opengeosys.org/df/d35/group__ODESolver.html).
 - Move up common Process parts from particular GroundwaterFlow process
   implementation. #951, #982
 - Separate Dirichlet boundary condition class implementation. #963
 - Split process output and post timestep. #998
 - Added pre- and postTimestep and -Iteration hooks to processes. #1094, #1100, #1101
 - New configuration tree parser
   - Checks configuration parameters more strictly, automatically prints error/warning messages.
   - Requires Boost >= 1.56 because of boost::optional with move semantics.
   - Command line argument `--config-warnings-nonfatal` that keeps OGS from terminating on warnings during
     configuration file parsing (errors still makes it terminate).
 - Axis aligned bounding box:
   - Is now a from the right half-open interval.
   - Removed template from class declaration.
 - MeshLib: Class MeshElementGrid implements a grid data structure supporting search operations.
 - Added cmake option `OGS_EIGEN_DYNAMIC_SHAPE_MATRICES` that makes OGS use dynamically.
   allocated shape matrices.
 - Added several cmake options for selecting which element types, dimensions and
   orders to be built. Selecting only few element types speeds up compilation
   significantly. #1092
 - Command line argument `-l` for OGS cli and testrunner to specify verbosity of logging. #1056
 - Added possibility to specify after which timesteps there shiuld be output.
 - Added possibility to specify timesteps of different size for use with
   transient processes.

#### DataExplorer and utilities
 - Added command line tool for creating layered meshes from raster files
 - OGSFileConverter is now a separate library
 - Raster file to structured grid conversion can now convert pixel values in user-defined scalar arrays
 - All scalar arrays will be displayed in mesh information window in DataExplorer
 - Added generation of structured meshes to DataExplorer
 - Restructured mesh creation access in DataExplorer
 - Mesh layers can be added to existing meshes in DataExplorer
 - Rework tools:
   - CreateBoundaryConditionsAlongPolyline
   - AddTopLayer
   - ResetPropertyInPolygonalRegion
   - removeMeshElements

### Infrastructure

- Minimum Boost version: 1.56.0. #943
- Boost requirement is now header-only. #940
- Optional support for VTK 7. #1083
- Test data is now a git submodule. #1000
- In-code defined Jenkins jobs. #970
- Use [clang's address and undefined behaviour sanitizers](https://svn.ufz.de:8443/job/OGS-6/job/Docker/job/clang-sanitizer/) on Jenkins now. #958


### Documentation

- Speed up builds with [ccache](http://docs.opengeosys.org/docs/devguide/advanced/using-ccache), #938
- Overview of the new non-linear, transient solver in [ODESolver](see
  http://doxygen.opengeosys.org/df/d35/group__ODESolver.html) source code
  documentation.

### Fixes
 - Fix bugs in GeoLib:
   - lineSegmentIntersects.
   - Polygon::splitPolygonAtIntersection.
   - Grid.
 - GeoMapper: Refactoring few methods, c++11. #977
 - Rework FileIO::GMSH interface
   - Process geometries located other than in the x-y-plane.
   - Respect the scaling factor for Stations.
   - Fix memory leaks.
   - Added/modified tests for GML-, GMS- and TetGen-files.

# 6.0.4

### Features:
 - Parallel computing framework for FEM by using PETSc, which also includes
   - Parallel input of partitioned mesh data.
   - Parallel output of solutions by using pvtu data format.
 - New data structures for mesh properties are used everywhere replacing
   Element's value member.
 - The penalty method to impose first-type boundary conditions was substituted
   with a non-penalty method for LIS and Eigen linear solvers.
 - Support for multiple nodal variables is extended to the boundary conditions,
   the sparsity pattern.
 - Passing of linear solver options from the project files is now possible.
 - The global matrix and global vector type of indices is consistent with the
   linear solver library being used.

### Infrastructure

- Added CMake option `OGS_CPU_ARCHITECTURE`, #858, [downloadable binaries](http://docs.opengeosys.org/download) build by Jenkins should now run on more CPUs
- Added CMake options for Boost, VTK and Eigen (`OGS_LIB_BOOST`, ...) to specify if libs are searched on the system first, then build them locally (`Default`), or you can specify to just use system libs (`System`) or force a local build (`Local`)
- Added CMake options for enabling Clang sanitizer:
  - `OGS_ADDRESS_SANITIZER`
  - `OGS_UNDEFINED_BEHAVIOR_SANITIZER`
- The zlib library is removed from ThirdParty directory.
- A LIS solver interface using Eigen's sparse matrices is now available through
  CMake option `OGS_USE_EIGENLIS`.
- CMake configuration uses [ccache](https://ccache.samba.org/) if available.

### Documentation

- Added [offline viewable Doxygen documentation](http://docs.opengeosys.org/docs/devguide/documentation/offline-documentation-viewer)

### Fixes
 - Fix all ogs-internal warnings on all OS.
 - Move eigen solver compute call to solve(); different fix for 0237275

# 6.0.3

### Features:
 - Mesh properties are now used for:
   - heterogeneous "initial conditions" (actually a start solution vector for the elliptic problem).
   - spatially heterogeneous hydraulic conductivity values in the groundwater flow process.
 - First steps towards time dependent problems: Time loop integration for processes is provided.
 - Interpolation of nodal quantities on elements using shape functions.

- Mesh generator can create surface meshes according to a given function

 - Utilities:
  - MoveGeometry

 - The DOF table handles now all of the provided element types: Hex 8 and 20, Line 2 and 3, Prism 6 and 15, Pyramid 5 and 13, Quad 4, 8, and 9, Tet 4 and 10, Triangle 3 and 6.
 - Eigen linear solver library can be used for solution of the linear systems of equations.

- Implemented [OctTree](https://github.com/ufz/ogs/pull/714) for fast searching
  points and nodes
 - Volumetric and surface grid
 - ElementSearcher NodeSearcher improvements
- Generalized the computation of rotation matrix to xy

### Fixes
 - FEFLOW interface supports element sets now.
 - Reduce compilation times by using forward declarations and removing unnecessary includes and using explicit template instantiation for often required classes.
 - GMSH2OGS: fixed bug in cases GMSH mesh does not contain line elements
 - CreateBoundaryConditionsAlongPolylines: fixed bug concerning the GeoLib and point ids.
 - PointVec corrected point id map
 - Shape interface creates polylines in a consistent state


### Infrastructure
- Replace quickcheck with autocheck. See https://github.com/thejohnfreeman/autocheck.git
for more details on autocheck
- Added support for cross-compiling with [MXE](http://mxe.cc/): build native Windows binaries on Linux and Mac OS, see [Cross-Compiling help page](http://docs.opengeosys.org/docs/devguide/advanced/cross-compiling) and #767
- Migrated to new Travis infrastructure (faster build times), see #775
- Simplified CMake library linking, see #769


## Test examples
- Test case: groundwater flow in the Unstrut catchment (model consists approximately of 9e6 hexahedral cells)
  - Simulations using homogeneous and heterogeneous hydraulic conductivity
  - Integrated rivers as Dirichlet type boundary conditions
  - Integrated groundwater recharge (spatialy homogeneous) as Neumann boundary condition

![unstrut_heterogeneous_rivers_top_layer_diff_recharge-no_recharge](https://cloud.githubusercontent.com/assets/329493/10043383/f75b837a-61f2-11e5-952e-f75d5d5a195b.png)

## Next steps
 The next big step will be the implementation of a parallelization scheme using PETSc library

### In development
- OGS#PETSc interface for parallel computing
- Solving of time dependent problems

### Planned
- Implementation of a linear parabolic pde solver
- Extending the linear elliptic solver to non-linear problems

# 6.0.2

| Released on 2015/06/15, [GitHub Release Link](https://github.com/ufz/ogs/releases/tag/6.0.2)

## Release notes

The second release of ogs6 introduces Neumann boundary conditions and VTK result output.

### Features:

- [Neumann boundary conditions](http://docs.opengeosys.org/docs/benchmarks/elliptic/groundwater-flow-neumann)
- Implement mesh properties for storage of data fields. This also includes mapped values (e.g. based on material id) [PR #542](https://github.com/ufz/ogs/pull/542), [PR #644](https://github.com/ufz/ogs/pull/644),
- Refactored mesh property classes to [enable VTK output](https://github.com/ufz/ogs/pull/692)
- Extended the available elements to quadratic (e.g. Quad8) based on generalized "element rules".  [#572](https://github.com/ufz/ogs/pull/572), [PR #656](https://github.com/ufz/ogs/pull/656), [PR #657](https://github.com/ufz/ogs/pull/657),
- Extend computation of shape matrices to lower dimensional elements embedded in higher dimensional space [#655](https://github.com/ufz/ogs/pull/655)
- Builds with MinGW (GCC) on Windows, see [Developer Guide](http://docs.opengeosys.org/docs/devguide/getting-started/prerequisites) and the new MinGW platform instructions
- [Cross-compiling for Windows with MXE on Mac OS](http://docs.opengeosys.org/docs/devguide/advanced/cross-compiling)
- [Support for new cross-platform IDE CLion](http://docs.opengeosys.org/docs/devguide/development-workflows/development-ides#clion)
- Add gradual refinement to the structured mesh generator [PR #539](https://github.com/ufz/ogs/pull/539)
- Add a command line tool "queryMesh" to search mesh information [PR #665](https://github.com/ufz/ogs/pull/665)
- Add a command line tool "AddTopLayer" to add an additional top layer (for example a soil layer, see also the [documentation](http://docs.opengeosys.org/docs/tools/meshing/addtoplayer)) [PR #649](https://github.com/ufz/ogs/pull/649)

### Fixes

- Performance optimizations in VTK mesh conversion, [PR #695](https://github.com/ufz/ogs/pull/695)
- Improve layered prism mesh construction
- Fix a lot of warnings on gcc/clang/msvc compilers improving the code


### Infrastructure

- Test runtime is monitored at Jenkins for [normal](https://svn.ufz.de:8443/job/OGS-6/job/Tests_Linux/) and [nightly large tests](https://svn.ufz.de:8443/job/OGS-6/job/Tests_Linux_Large/)
- Utilities are build by separate Jenkins Jobs, e.g. [Win_Utils](https://svn.ufz.de:8443/job/OGS-6/job/Win_Utils/)

## Test examples

- [Groundwater flow (Neumann)](http://docs.opengeosys.org/docs/benchmarks/elliptic/groundwater-flow-neumann):
![](http://docs.opengeosys.org/assets/files/Documentation/Selected-Benchmarks/groundwaterflow-neumann/square_1e2_neumann_result.png)

## Next steps

### In development

- Heterogeneous fields (for e.g. hydraulic conductivity parameters)
- Octree data structures for fast searches
- OGS#PETSc interface for parallel computing

### Planned

- Parallelization scheme using PETSc library
- Extending the linear elliptic solver to non-linear problems


--------------------------------------------------------------------------------

# 6.0.1

| Released 2015/03/02, [GitHub Release Link](https://github.com/ufz/ogs/releases/tag/6.0.1)

The 6th version of OpenGeoSys (OGS) is under way. After single and coupled FORTRAN modules in ROCKFLOW 1+2, the C version 3 with dynamic data structures, the object-oriented C++ parallelized version 4, completed with data integration and visualization tools by version 5; ogs6 - as an open source project - is aimed at performing on supercomputing platforms and providing complete workflows for solving of coupled multi-field problems in real world applications. The major paradigms of ogs6 are being developer-friendly, performing, and user-friendly.

## Important links:

- Getting started tutorial: http://docs.opengeosys.org/docs/quickstart
- Descriptions of selected benchmarks: http://docs.opengeosys.org/docs/benchmarks

- Source code access: https://github.com/ufz/ogs
- Developer guide: http://docs.opengeosys.org/docs/devguide

## Release notes
The first version ogs6 is dedicated for elliptic problems.

### Features:
 - Basic structures of processes
 - Mathematical operations are based on Eigen3 library
 - Linear solvers: DenseMatrix with Gauss elimination, and LIS (http://www.ssisc.org/lis/)
 - XML based IO
 - Standard finite element method (FEM)
 - Available element types: lines, triangles, quads, hexahedra
 - Dirichlet boundary conditions
 - Linear elliptic solver (e.g. Groundwater flow) for scalar quantities in homogeneous media

### Fixes
- DenseMatrix Gauss algorithm pivoting
- Fixing mem-leaks on DataExplorer start up
- Fixing resizing and layout issues in various DataExplorer dialogs

## Test examples
![](https://cloud.githubusercontent.com/assets/329493/6170573/ce9fd96c-b2d5-11e4-9936-a470e7be281f.png)

 - Example 1: Unit square (access)
 - Example 2: Unit cube (access): has 1000 hexahedra  elements with Dirichlet boundary conditions (u=1|x=0) and (u=-1|x=1)

## Next steps

### In development
 - OGS#PETSc interface for parallel computing (02/2015*planned)
 - Neuman boundary conditions (03/2015*planned)

### Planned
 - Parabolic solver
