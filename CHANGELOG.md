# 6.2.0

## Features

### New processes

- HeatTransportBHE process supporting 1U, CXA, and CXC BHE types. [#2221](https://github.com/ufz/ogs/pull/2221), [#2332](https://github.com/ufz/ogs/pull/2332),
  [#2271](https://github.com/ufz/ogs/pull/2271), [#2275](https://github.com/ufz/ogs/pull/2275)
- Staggered implementation of a thermo-mechanical with phasefield process. [#2102](https://github.com/ufz/ogs/pull/2102)
- Richards mechanics process. [#2189](https://github.com/ufz/ogs/pull/2189)
- Small deformation process with non-local integration of damage. [#2294](https://github.com/ufz/ogs/pull/2294)
- Staggered implementation of phasefield process. [#2052](https://github.com/ufz/ogs/pull/2052)
- ComponentTransport process in revised formulation. [#2200](https://github.com/ufz/ogs/pull/2200)
- Multi-component transport process. [#2304](https://github.com/ufz/ogs/pull/2304)

### Other process' changes

- A Jacobian tester: a process's Jacobian assembly can be compared to a
  numerical Jacobian (mostly for development. [#2238](https://github.com/ufz/ogs/pull/2238)
- Add `setInitialConditions()` call to processes and local assemblers. [#2334](https://github.com/ufz/ogs/pull/2334)
- Several bug fixes for LIE/HM process including "fracture into matrix
  leak-off", Darcy velocity output in the fracture. [#2129](https://github.com/ufz/ogs/pull/2129)
- Support for intersecting fractures (x-crossing and t-junction) in LIE/SD
  process. [#2235](https://github.com/ufz/ogs/pull/2235), [#2293](https://github.com/ufz/ogs/pull/2293)
- Fixed the calculation of the Darcy velocity in staggered TH. [#2127](https://github.com/ufz/ogs/pull/2127)

#### Numerics

- Staggered scheme for coupled processes with different orders of elements.
  [#2016](https://github.com/ufz/ogs/pull/2016)
- Subdomain deactivation within time intervals. [#2297](https://github.com/ufz/ogs/pull/2297)
- Add a driver for an iteration based time stepping algorithm. [#2318](https://github.com/ufz/ogs/pull/2318)

#### Boundary condition

- Implementation of Python boundary conditions. [#2170](https://github.com/ufz/ogs/pull/2170)
- Implementation of constraint boundary conditions. [#2145](https://github.com/ufz/ogs/pull/2145)
- Dirichlet boundary condition within a time interval. [#2272](https://github.com/ufz/ogs/pull/2272)
- BoundaryElementSearch: Return bulk element id and bulk element face id. [#2125](https://github.com/ufz/ogs/pull/2125)
- Removed Neumann boundary condition for displacement jumps in LIE processes.
  [#2153](https://github.com/ufz/ogs/pull/2153)

#### Source term

- Use parameter for source terms. [#2061](https://github.com/ufz/ogs/pull/2061)
- Volumetric source terms implementation. [#2220](https://github.com/ufz/ogs/pull/2220), [#2234](https://github.com/ufz/ogs/pull/2234), [#2241](https://github.com/ufz/ogs/pull/2241), [#2261](https://github.com/ufz/ogs/pull/2261)

#### Input and output

- Writing and reading of integration point data. Implemented sigma and epsilon
  output for some processes. [#2071](https://github.com/ufz/ogs/pull/2071), [#2203](https://github.com/ufz/ogs/pull/2203), [#2324](https://github.com/ufz/ogs/pull/2324)
- Add input of vtu-meshes for boundary conditions additionally to the gml input.
  This is later used by the heterogeneous parameters and source terms. [#2140](https://github.com/ufz/ogs/pull/2140),
  [#2141](https://github.com/ufz/ogs/pull/2141), [#2156](https://github.com/ufz/ogs/pull/2156)
- Parameters may now be explicitly defined on arbitrary subdomains. This merges
  the Heterogeneous Dirichlet and Neumann boundary conditions with their,
  previously only homogeneous, counterparts. [#2376](https://github.com/ufz/ogs/pull/2376)
- Parameters now support space-dependent function input via exprtk library.
  [#2325](https://github.com/ufz/ogs/pull/2325), [#2339](https://github.com/ufz/ogs/pull/2339)
- Output of primary variables on arbitrary subdomains. [#2372](https://github.com/ufz/ogs/pull/2372), [#2299](https://github.com/ufz/ogs/pull/2299)
- Output is possible at specific times for adaptive time stepping and
  evolutionaryPIDController. [#2079](https://github.com/ufz/ogs/pull/2079)
- Calculate and output specific flux. [#2411](https://github.com/ufz/ogs/pull/2411)
- Enable surface flux calculation for component transport process. [#2168](https://github.com/ufz/ogs/pull/2168)
- Interpolated pressure (on higher order elements' nodes) output for
  Richards-mechanics and hydro-mechanics processes. [#2228](https://github.com/ufz/ogs/pull/2228)
- Improve output of nodal aperture and aperture vector in LIE/HM. [#2050](https://github.com/ufz/ogs/pull/2050)
- Add output of nodal forces and hydraulic flow in mechanics and coupled
  mechanics processes, SD, HM, LIE/SD, LIE/HM. [#2118](https://github.com/ufz/ogs/pull/2118)
- Enable surface flux calculation for HT process. [#2132](https://github.com/ufz/ogs/pull/2132)
- Fixed pvd output. [#2036](https://github.com/ufz/ogs/pull/2036)

### Material models

- BGRa creep model. [#2167](https://github.com/ufz/ogs/pull/2167)
- New cohesive zone mode I fracture model for LIE processes. [#2142](https://github.com/ufz/ogs/pull/2142), [#2157](https://github.com/ufz/ogs/pull/2157)
- Add MFront/TFEL solid constitutive relation support via.
  MFrontGenericInterfaceSupport library. CMake option `OGS_USE_MFRONT`. [#2259](https://github.com/ufz/ogs/pull/2259)
- Infrastructure for multi-phase, multi-component material properties library.
  [#2303](https://github.com/ufz/ogs/pull/2303),
  [#2415](https://github.com/ufz/ogs/pull/2415)
- Anisotropic tensors may now be given in given local coordinate system. [#2370](https://github.com/ufz/ogs/pull/2370)
- Non-constant density model implementation in HC process. [#2200](https://github.com/ufz/ogs/pull/2200)
- Add second derivatives of permeability functions in Richards flow. [#2188](https://github.com/ufz/ogs/pull/2188)
- Different solid material models can now be defined on different materialIDs.
  [#2216](https://github.com/ufz/ogs/pull/2216), [#2262](https://github.com/ufz/ogs/pull/2262), [#2270](https://github.com/ufz/ogs/pull/2270)
- Move solid constitutive relation creation in single place. [#2160](https://github.com/ufz/ogs/pull/2160)

### Testing and documentation

- Migrated Appveyor tests to [Azure
  Pipelines](https://dev.azure.com/ogsci/ogs/_build). [#2342](https://github.com/ufz/ogs/pull/2342)
- Added cppcheck, clang-tidy and include-what-you-use. [#2078](https://github.com/ufz/ogs/pull/2078), [#2328](https://github.com/ufz/ogs/pull/2328), [#2377](https://github.com/ufz/ogs/pull/2377)
- Added check for header standalone compilation, can be enabled with
  `OGS_CHECK_HEADER_COMPILATION=ON`. [#2043](https://github.com/ufz/ogs/pull/2043)
- Jenkins shows nice summaries of compiler warnings. [#2206](https://github.com/ufz/ogs/pull/2206)
- Large tests are fixed and run on Jenkins upon master merge. [#2056](https://github.com/ufz/ogs/pull/2056), [#2155](https://github.com/ufz/ogs/pull/2155)
- Re-enabled code coverage reports (for the testrunner only) with
  [Codecov](http://codecov.io/gh/ufz/ogs). [#2333](https://github.com/ufz/ogs/pull/2333), [#2336](https://github.com/ufz/ogs/pull/2336)
- Commits containing `[ci skip]` in the commit message do not trigger a Jenkins
  build.
- Add `vtkdiff` test configuration to project files, s.t. the vtkdiff tests are
  performed after successful run comparing output to reference files. This
  possibility is also reflected in a new CMake function `OgsTest` as an
  alternative to the `AddTest`. [#2255](https://github.com/ufz/ogs/pull/2255), [#2257](https://github.com/ufz/ogs/pull/2257)

### New tools

- `TecPlotTools`: splitting files containing multiple zones into seperate
  TecPlot files. [#2114](https://github.com/ufz/ogs/pull/2114)
- `TecPlot-Reader`: converting TecPlot rasters into OGS meshes (one file per
  zone, containing all variables as scalar arrays). [#2114](https://github.com/ufz/ogs/pull/2114)
- [`constructMeshesFromGeometry`](https://www.opengeosys.org/docs/tools/model-preparation/constructmeshesfromgeometry/):
  Construction of boundary meshes from bulk mesh and gml files. [#2144](https://github.com/ufz/ogs/pull/2144)
- [`identifySubdomains`](https://www.opengeosys.org/docs/tools/model-preparation/identifysubdomains/):
  Identification of boundary meshes (or any subdomains in general) in the bulk
  mesh. Performs geometrical tests and creates and verifies necessary
  `bulk_node_ids` and `bulk_element_ids` maps. [#2227](https://github.com/ufz/ogs/pull/2227), [#2252](https://github.com/ufz/ogs/pull/2252)
- `Mesh2Raster`: converts 2D OGS meshes into raster files of arbitrary pixel
  size, where node elevation is represented by pixel value. [#2367](https://github.com/ufz/ogs/pull/2367)
- `GocadSGridReader` tool reading the Gocad/SKUA stratigraphic grid format and
  writing the data in the vtu format. [#2316](https://github.com/ufz/ogs/pull/2316)

### New tools and CLI usage

- `ogs --help` shows the given CMake options. [#2210](https://github.com/ufz/ogs/pull/2210)
- Unify command line version info output. [#2194](https://github.com/ufz/ogs/pull/2194)
- Rewrite `partmesh` tool and add partitioning of boundary meshes (or subdomains
  in general) according to the partition of the bulk mesh. [#2159](https://github.com/ufz/ogs/pull/2159), [#2178](https://github.com/ufz/ogs/pull/2178), [#2195](https://github.com/ufz/ogs/pull/2195)
- Add new features to `ExtractSurface` tool. [#2387](https://github.com/ufz/ogs/pull/2387), [#2401](https://github.com/ufz/ogs/pull/2401)
- updated utility `moveMeshNodes`: algorithm for mesh on mesh mapping now
  calculates exact node elevation instead of using interpolation. [#2390](https://github.com/ufz/ogs/pull/2390)

### Data Explorer

- Listing of source terms and boundary conditions in Data Explorer DataView (no
  visualisation yet). [#2110](https://github.com/ufz/ogs/pull/2110)
- Mesh element removal can now remove elements based on value ranges of
  arbitrary scalar arrays (currently only int- and double arrays are supported).
  [#2115](https://github.com/ufz/ogs/pull/2115)
- added custom VTK filter to represent raster data as point clouds. [#2121](https://github.com/ufz/ogs/pull/2121)
- geometrical points can now be converted into station points. [#2369](https://github.com/ufz/ogs/pull/2369)
- fixed issue with geometrical surfaces not being loaded correctly. [#2388](https://github.com/ufz/ogs/pull/2388)
- Replace deprecated QVTKWidget with QVTKOpenGLWidget [#2432](https://github.com/ufz/ogs/pull/2432)

### Other notable code changes

- C++17 standard is enabled and is allowed in the production code (given the
  compiler support). [#2298](https://github.com/ufz/ogs/pull/2298)
- Separate monolithic ProcessLib into individual processes. Now it is possible
  to build ogs with selected processes only. This also improves linking times.
  [#2017](https://github.com/ufz/ogs/pull/2017)
- Parameters are extracted in own library. [#2413](https://github.com/ufz/ogs/pull/2413)
- Port secondary variable extrapolation and output for PETSc builds. [#2082](https://github.com/ufz/ogs/pull/2082)
- Extend Kelvin mapping functions and move implementation to MathLib. [#2060](https://github.com/ufz/ogs/pull/2060),
  [#2075](https://github.com/ufz/ogs/pull/2075), [#2044](https://github.com/ufz/ogs/pull/2044)
- Collect generic algorithms in single header file. [#2161](https://github.com/ufz/ogs/pull/2161)
- Remove unused MeshSubsets class. [#2135](https://github.com/ufz/ogs/pull/2135)
- Removed writing of xsd header in XML files, [#2198](https://github.com/ufz/ogs/pull/2198)

## Infrastructure

- Migrated LFS storage from GitLab to [Artifactory](https://ogs.jfrog.io/ogs).
  [#2359](https://github.com/ufz/ogs/pull/2359)
- Optimized ctest runtime by starting long-running benchmarks first. [#2310](https://github.com/ufz/ogs/pull/2310)
- Proper RPATH handling for shared library installations. [#2208](https://github.com/ufz/ogs/pull/2208)
- [Package OGS inside
  container](https://www.opengeosys.org/docs/userguide/basics/container/) with
  [Singularity](https://www.sylabs.io/singularity/); [more
  docs](https://www.opengeosys.org/docs/devguide/advanced/singularity/). [#2193](https://github.com/ufz/ogs/pull/2193),
  [#2356](https://github.com/ufz/ogs/pull/2356)
- Migrated opengeosys.org to a static site generator ([Hugo](https://gohugo.io))
  unifying documentation and general OGS info. [#2088](https://github.com/ufz/ogs/pull/2088), [#2095](https://github.com/ufz/ogs/pull/2095), [#2123](https://github.com/ufz/ogs/pull/2123)
- Speed-up CMake run time. [#2072](https://github.com/ufz/ogs/pull/2072), [#2392](https://github.com/ufz/ogs/pull/2392)

### CMake options changes

- `OGS_USE_PYTHON` enables Python BCs. [#2170](https://github.com/ufz/ogs/pull/2170)
- `OGS_BUILD_TESTS` was renamed to `BUILD_TESTING`. [#2350](https://github.com/ufz/ogs/pull/2350)
- Added `OGS_USE_CVODE`. [#2344](https://github.com/ufz/ogs/pull/2344)
- Added `OGS_BUILD_PROCESSES` for `;`-separated list of processes to build.
  [#2233](https://github.com/ufz/ogs/pull/2233)
- `OGS_USE_CONAN=ON` is now the default when `conan` was found. [#2207](https://github.com/ufz/ogs/pull/2207)

### Version info

- CMake minimum version 3.10
- Visual Studio minimum (and tested) version 2017
- GCC minimum version 6.2 (tested: 6.4)
- Clang minimum version 3.5 (tested: 7.0)
- Boost minimum version 1.66.0
- VTK minimum version 8.1. [#2158](https://github.com/ufz/ogs/pull/2158)
- Qt tested version 5.11.2
- Python tested version 3.7.2


# 6.1.0

The changes since the prerelease 6.1.0-rc1 contain few bug fixes and
improvements of the CI.

### Features:

#### New processes:
- ComponentTransport https://github.com/ufz/ogs/pull/1758
- PhaseField https://github.com/ufz/ogs/pull/1813 and https://github.com/ufz/ogs/pull/1895
- RichardsComponentTransport https://github.com/ufz/ogs/pull/1929
- ThermoMechanics https://github.com/ufz/ogs/pull/1800
- TwoPhaseFlow p-rho https://github.com/ufz/ogs/pull/1613

#### Other process' changes:
- New equation assembly approach for the staggered scheme. With that,  the
  coupling assembly computations are  changed from  performing assembly across
  the different Process classes for a coupling to a single Process class for
  coupled processes, which was assumed only for the monolithic before. With the
  changes, the original Process class that for monolithic scheme originally can
  now handle both of the monolithic and staggered schemes.  So far,  HT classes
  get the staggered scheme based on this new framework of  assembly.
  https://github.com/ufz/ogs/pull/1970
- Heterogeneous liquid flow properties  (https://github.com/ufz/ogs/pull/1979,
  https://github.com/ufz/ogs/pull/2017)
- New boundary conditions added: Nonuniform Dirichlet
  (https://github.com/ufz/ogs/pull/1952) and Neumann
  (https://github.com/ufz/ogs/pull/1891), and normal traction
  (https://github.com/ufz/ogs/pull/1896)
- Framework for time stepping and a first application of adaptive time stepping,
  EvolutionaryPIDcontroller, and automatic time step control
  (https://github.com/ufz/ogs/pull/1803).
- Nodal source terms (https://github.com/ufz/ogs/pull/1977)
- Fix deformation processes to work correctly with PETSc parallelization
  (https://github.com/ufz/ogs/pull/1838).
- Fix the access to local data of PETScVector (https://github.com/ufz/ogs/pull/1797).
- Add damping factor to global Newton. https://github.com/ufz/ogs/pull/2004
- Extend extrapolator to vectorial quantities and replace the component-wise
  output of Darcy velocity and the stress/strain in mechanical processes with
  single vectorial output.

#### Material model changes
- Separate FractureModels in LIE https://github.com/ufz/ogs/pull/1971
- Add material forces as published in
  http://www.sciencedirect.com/science/article/pii/S0093641317303865
  https://github.com/ufz/ogs/pull/1936

#### Testing and documentation:
- Insitu visualization with ParaView Catalyst. See
  [presentation](https://github.com/ufz/ogs/files/867280/Insitu-Department.pdf).
  [#1744](https://github.com/ufz/ogs/pull/1744), [#1732](https://github.com/ufz/ogs/pull/1732). As a consequence VTK 7.1 is now required.
- Benchmark docs are now part of the code (in `web/content`) and can contain
  [interactive 3D
  visualizations](https://dev.opengeosys.org/docs/benchmarks/elliptic/groundwater-flow-neumann/#results-and-evaluation)
  via [vtk.js](https://kitware.github.io/vtk-js/). [#1706](https://github.com/ufz/ogs/pull/1706), [#1714](https://github.com/ufz/ogs/pull/1714), [#1723](https://github.com/ufz/ogs/pull/1723), [#1729](https://github.com/ufz/ogs/pull/1729).
- Migrated handling of test data files from *git-submodule* to *git-lfs*, see
  [docs](https://docs.opengeosys.org/docs/devguide/testing/test-data). [#1964](https://github.com/ufz/ogs/pull/1964),
  [#1982](https://github.com/ufz/ogs/pull/1982), [#1984](https://github.com/ufz/ogs/pull/1984), [#2010](https://github.com/ufz/ogs/pull/2010), [#2012](https://github.com/ufz/ogs/pull/2012).  Now  [git-lfs](https://git-lfs.github.com/) is
  required. Check the
  [installation](https://docs.opengeosys.org/docs/devguide/getting-started/prerequisites)
  instructions.

### Infrastructure:
- Fully moved to Conan for automatic third-party library handling. Can be
  enabled with `OGS_USE_CONAN=ON`, see
  [docs](https://docs.opengeosys.org/docs/devguide/advanced/conan-package-manager).
  [#1907](https://github.com/ufz/ogs/pull/1907)
- Conan version 1.0 is now required.
- Dropped travis CI environment and added few new tests on Jenkins because of
  simpler maintenance.

#### CMake options changes:
- `OGS_EIGEN_DYNAMIC_SHAPE_MATRICES` defaults to OFF on Release
  config, ON otherwise. Can be overridden by explicitly setting the option. [#1673](https://github.com/ufz/ogs/pull/1673)
- New `OGS_EIGEN_INITIALIZE_MATRICES_BY_NAN` defaults to ON for easier spotting
  of non-initialized matrices. When OFF, the Eigen's default initialization to 0
  is skipped resulting in slightly faster execution.
  https://github.com/ufz/ogs/pull/1917
- Set default Eigen's cmake flag disabling vectorization since this lead to
  several problems in different environments.
  https://github.com/ufz/ogs/pull/1919 and the issue linked there
  https://github.com/ufz/ogs/issues/1881

#### Other
- PETSc config is tested on Jenkins (envinf1)
- OGS binaries are provided as eve / envinf1 modules. See
  [docs](https://docs.opengeosys.org/docs/quickstart/basics/envinf1) for
  details. [#1753](https://github.com/ufz/ogs/pull/1753)
- Migrated Data Explorer to Qt5. [#1622](https://github.com/ufz/ogs/pull/1622), [#1625](https://github.com/ufz/ogs/pull/1625)
- Windows builds are tested on MSVC 2017 on own hardware and on MS


# 6.0.8

The highlight of the release is the implementation of the Lower-Interface
Elements for both the small deformation process (M) and hydro-mechanics process
(HM) allowing fractures to be incorporated in the solution domain.
For the liquid flow and two-phase flow processes several material models for the
fluids pressure, density, permeability, and viscosity were added.

### Features:

 - Implementation of hydro-mechanics (HM) with LIE. [#1537](https://github.com/ufz/ogs/pull/1537)-[#1541](https://github.com/ufz/ogs/pull/1541)
 - Implementation of small deformation (M) with LIE. [#1452](https://github.com/ufz/ogs/pull/1452)
 - Fracture constitutive models. [#1434](https://github.com/ufz/ogs/pull/1434)
 - Hydro-Mechanics process. [#1508](https://github.com/ufz/ogs/pull/1508)
 - First version of monolithic hydro-thermal process implementation with
   Boussinesq approximation using constant viscosity. [#1534](https://github.com/ufz/ogs/pull/1534)
 - Two phase flow process with pp model implementation. [#1530](https://github.com/ufz/ogs/pull/1530)
 - Richards flow process implementation. [#1473](https://github.com/ufz/ogs/pull/1473)
 - Liquid process. [#1468](https://github.com/ufz/ogs/pull/1468)
 - Classes for relative permeability models. [#1531](https://github.com/ufz/ogs/pull/1531)
 - Classes for capillary models. [#1517](https://github.com/ufz/ogs/pull/1517), [#1578](https://github.com/ufz/ogs/pull/1578)
 - Ehlers single-surface yield function constitutive relation model. [#1556](https://github.com/ufz/ogs/pull/1556)
 - Support scaling, GMRES, and Pardiso in Eigen linear solvers. [#1509](https://github.com/ufz/ogs/pull/1509) [#1510](https://github.com/ufz/ogs/pull/1510)
 - Piecewise linear Monotonic curve and a generic curve parser. [#1529](https://github.com/ufz/ogs/pull/1529)
 - Support searching boundary nodes in MeshLib::NodeSearch. [#1459](https://github.com/ufz/ogs/pull/1459)
 - Support specifying the shape function order in process variables. [#1503](https://github.com/ufz/ogs/pull/1503)
 - Command line option --unbuffered-std-out to deactivate buffer for standard output. [#1514](https://github.com/ufz/ogs/pull/1514)
 - CMake option OGS_FATAL_ABORT for debugging. [#1432](https://github.com/ufz/ogs/pull/1432)
 - Set the default OGS_LOG_LEVEL to debug in release builds. [#1522](https://github.com/ufz/ogs/pull/1522)
 - Add integration order in input files. [#1464](https://github.com/ufz/ogs/pull/1464)
 - Migrated Data Explorer to Qt5. [#1622](https://github.com/ufz/ogs/pull/1622), [#1625](https://github.com/ufz/ogs/pull/1625)
 - Benchmarks can be run on specific configurations only by using the new parameter
   `REQUIREMENTS` in `AddTest()` (in CMake). [#1610](https://github.com/ufz/ogs/pull/1610)

### Utilities
New utilities:
 - createQuadraticMesh [#1500](https://github.com/ufz/ogs/pull/1500)
 - convertToLinearMesh [#1554](https://github.com/ufz/ogs/pull/1554)
 - postLIE [#1555](https://github.com/ufz/ogs/pull/1555)

New features:
 - extend NodeReordering to correct ordering of nonlinear nodes [#1519](https://github.com/ufz/ogs/pull/1519)


### Infrastructure:

 - Ctest now works on Windows too by removing time-wrappers. [#1480](https://github.com/ufz/ogs/pull/1480)
 - Moved to public Jenkins instance at jenkins.opengeosys.org. [#1505](https://github.com/ufz/ogs/pull/1505)
 - Doxygen warnings parser in Jenkins will mark a build as unstable
   if there are Doxygen warnings. [#1585](https://github.com/ufz/ogs/pull/1585)
 - Benchmarking on Jenkins now saves the standard output into a file for each
   test. [#1528](https://github.com/ufz/ogs/pull/1528)

### Fixes:
 - Fix LocalToGlobalIndexMap with mutliple variables and with multiple componets. [#1433](https://github.com/ufz/ogs/pull/1433) [#1440](https://github.com/ufz/ogs/pull/1440)
 - Fix PropertyVector<T*> for multi-component case. [#1441](https://github.com/ufz/ogs/pull/1441)
 - Fix checking IDs of nonlinear nodes. [#1495](https://github.com/ufz/ogs/pull/1495)
 - Fix incorrect use of getNumberOfBaseNodes(). [#1515](https://github.com/ufz/ogs/pull/1515)
 - Fix computing sparsity pattern for mixed shape function order cases. [#1548](https://github.com/ufz/ogs/pull/1548)
 - Fix that iterations and residuals were not printed when Eigen linear solver fails. [#1499](https://github.com/ufz/ogs/pull/1499)
 - Fix all of the Doxygen warnings in the code. [#1569](https://github.com/ufz/ogs/pull/1569) [#1573](https://github.com/ufz/ogs/pull/1573)
 - Fix all of the input file/keyword documentation and its generation.


# 6.0.7

### Features:

The main features of this release is the implementation of two new processes,
the small deformation, and the heat conduction.  Some extensions were done to
the DOF table to be able to manage multi-component/multi-variable processes.
Also, during implementation of the Robin boundary conditions, the base classes
for boundary conditions were generalized.

 - Add small deformation process with linear elastic material model. The
   implementation is based on the Kelvin mapping. [#1340](https://github.com/ufz/ogs/pull/1340)
 - Added B-Matrices and Kelvin mapping tools for deformation processes. [#1359](https://github.com/ufz/ogs/pull/1359)
 - Heat conduction process implementation. [#1328](https://github.com/ufz/ogs/pull/1328)
 - Finalize support for multicomponent boundary conditions adding configuration
   parser. [#1343](https://github.com/ufz/ogs/pull/1343)
 - Added uniform Robin boundary condition. [#1336](https://github.com/ufz/ogs/pull/1336)
 - Added a generic natural boundary condition class. [#1337](https://github.com/ufz/ogs/pull/1337)
 - Added Robin boundary condition. [#1336](https://github.com/ufz/ogs/pull/1336)
 - Reworked the Parameter class. It now serves as a basis for BCs and ICs.
   [#1357](https://github.com/ufz/ogs/pull/1357), [#1356](https://github.com/ufz/ogs/pull/1356)
 - Added time-dependent Dirichlet BCs. [#1380](https://github.com/ufz/ogs/pull/1380)
 - Add calculation of surface flux, tests for groundwater flow. [#1429](https://github.com/ufz/ogs/pull/1429)
 - Implemented numerical Jacobian assembly for Newton-Raphson solver. [#1352](https://github.com/ufz/ogs/pull/1352)
 - Added the new parameter type "Group" which can be used for setting material
   ID dependent values. [#1426](https://github.com/ufz/ogs/pull/1426)
 - Added fluid property class and several fluid density and viscosity models
   based on it. [#1398](https://github.com/ufz/ogs/pull/1398), [#1435](https://github.com/ufz/ogs/pull/1435)
 - Enabled solving of axially symmetric problems on 2D meshes for all currently
   implemented processes. [#1443](https://github.com/ufz/ogs/pull/1443)
 - Added time measurement for assembly and solvers. [#1322](https://github.com/ufz/ogs/pull/1322)
 - Added named functions, out of which expressions can be built up at run-time
   from the prj file, which can be used to output additional nodal quantities.
   [#1314](https://github.com/ufz/ogs/pull/1314), [#1315](https://github.com/ufz/ogs/pull/1315)
 - Added component-wise norms, and flexible convergence criteria for nonlinear
   solvers. [#1349](https://github.com/ufz/ogs/pull/1349), [#1342](https://github.com/ufz/ogs/pull/1342)
 - Restructured the time loop. [#1364](https://github.com/ufz/ogs/pull/1364)

### Utilities
New utilities:
 - createNeumannBc: The tool integrates the given element property and writes
   the computed data as a PropertyVector with the name `node_aggregated_gwn`
   into the mesh. The tool also outputs an OGS-5 direct source term (Neumann
   boundary condition) data file. [#1346](https://github.com/ufz/ogs/pull/1346)
 - scaleProperty for simple rescaling of mesh properties. [#1347](https://github.com/ufz/ogs/pull/1347)
 - convertGEO for geometric file conversion, e.g. gli to glm. [#1360](https://github.com/ufz/ogs/pull/1360)
 - swapNodeCoordinateAxes to swap node coordinate values, e.g. XY to XZ plane.
   [#1361](https://github.com/ufz/ogs/pull/1361)
New feature:
 - Support tetrahedra types in generateStructuredMesh. [#1353](https://github.com/ufz/ogs/pull/1353)

### Infrastructure:
 - Migrated all important Jenkins jobs to script-based
   [Jenkins Pipeline](https://jenkins.io/doc/pipeline/)
   functionality.  For an introduction see
   [docs.opengeosys.org - Continuous Integration](https://docs.opengeosys.org/docs/devguide/development-workflows/continuous-integration).
   [#1392](https://github.com/ufz/ogs/pull/1392), [#1396](https://github.com/ufz/ogs/pull/1396), [#1404](https://github.com/ufz/ogs/pull/1404), [#1411](https://github.com/ufz/ogs/pull/1411), [#1424](https://github.com/ufz/ogs/pull/1424), [#1428](https://github.com/ufz/ogs/pull/1428), [#1436](https://github.com/ufz/ogs/pull/1436)
 - Moved CMake logic for packaging executable dependencies (such as shared libs)
   to the install and package targets instead of running after each executable
   gets build. [#1458](https://github.com/ufz/ogs/pull/1458)
 - Increase minimum clang compiler version to 3.5 in course of updating travis
   build enviroment to Ubuntu LTS 14.04. [#1448](https://github.com/ufz/ogs/pull/1448)
 - Added a script that generates crosslinked Doxygen pages out of ctest input
   files [#1348](https://github.com/ufz/ogs/pull/1348)

### Fixes:
 - Fix an issue that a shape vector was defined as a column vector. Corrected to
   a row vector . [#1288](https://github.com/ufz/ogs/pull/1288)
 - Fix usage of `boost::optional<T const&>`, which has changed in version 1.61.
   [#1385](https://github.com/ufz/ogs/pull/1385)
 - Fix Grid (enlarge bounding box to fit all points). [#1369](https://github.com/ufz/ogs/pull/1369)
 - Fix mapping of geometries to mesh surfaces. [#1327](https://github.com/ufz/ogs/pull/1327). [#1368](https://github.com/ufz/ogs/pull/1368)
 - Fix transmitting raster data to element properties. [#1347](https://github.com/ufz/ogs/pull/1347)
 - Fixed missing XSD files in packages. [#1410](https://github.com/ufz/ogs/pull/1410)
 - Fix a shape vector to a row vector. [#1288](https://github.com/ufz/ogs/pull/1288)
 - Fix FEFLOW import. [#1397](https://github.com/ufz/ogs/pull/1397)
 - Fix NodeReordering to check ordering of each element. [#1425](https://github.com/ufz/ogs/pull/1425)
 - Fix builds linking shared VTK library. [#1431](https://github.com/ufz/ogs/pull/1431)
 - Fix global Newton iteration counter. [#1341](https://github.com/ufz/ogs/pull/1341)
 - Correct few loops over mesh nodes, which should run over the mesh subsets.
   [#1437](https://github.com/ufz/ogs/pull/1437)
 - Fix shape function computation for 2D elements lying in the x-y-plane [#1318](https://github.com/ufz/ogs/pull/1318)
 - Fix AddTest, s.t. ctest now really checks results. [#1325](https://github.com/ufz/ogs/pull/1325)
 - Made Eigen preconditioner configurable. [#1367](https://github.com/ufz/ogs/pull/1367)

# 6.0.6

### Features:
 - Add external ode-solver interface with [Sundials CVODE
   library](http://computation.llnl.gov/projects/sundials-suite-nonlinear-differential-algebraic-equation-solvers/cvode). [#1109](https://github.com/ufz/ogs/pull/1109)
 - Add piecewise linear curves parser to the project files. The curves are
   specified by two vectors, the coordinates and values. They can be used for
   example to map temporal dependencies (time-dependent boundary conditions) or
   as approximations of coefficient dependencies (e.g. pressure-saturation
   curves). [#1149](https://github.com/ufz/ogs/pull/1149)
 - Extend the LocalAssemblerInterface by adding default implementations of
   pre/postTimestep and assembleJacobian functions. The current time and time
   step size are passed in the preTimestep call to the particular processes. [#1214](https://github.com/ufz/ogs/pull/1214)
 - Add support multi-variable/multi-component in the DOF table interface and
   extend the initial conditions to multi-components. [#1224](https://github.com/ufz/ogs/pull/1224)
 - Major rework of the general process interface. That also affects the
   interface of concrete processes and local assemblers.
   least squares optimization. [#1145](https://github.com/ufz/ogs/pull/1145)
 - Added functionality for the output of secondary variables. [#1171](https://github.com/ufz/ogs/pull/1171)
 - Added material properties for zeolite adsorption and Calcium oxide/hydroxide
   reactions. [#1178](https://github.com/ufz/ogs/pull/1178)
 - Transferred the TES process, a monolithically coupled THC model for simulating
   thermochemical energy storag devices, from OGS5. [#1181](https://github.com/ufz/ogs/pull/1181)
 - Introduced a general scheme for documenting OGS6 input file settings. #978
 - Added copy constructor for the class Surface, minor improvements in GeoLib. [#1237](https://github.com/ufz/ogs/pull/1237)
 - Added classes GeoLib::LineSegment and GeoLib::Polyline::SegmentIterator. [#1139](https://github.com/ufz/ogs/pull/1139)
 - GMSHInterface can handle stations as constraints. [#1212](https://github.com/ufz/ogs/pull/1212)
 - Added functionality to duplicate geometric data. [#1229](https://github.com/ufz/ogs/pull/1229)
 - Station names can be modified in Data Explorer [#1273](https://github.com/ufz/ogs/pull/1273)

### Infrastructure
 - Fix circular dependencies on library level. This allows for dynamic linking
   which is faster than static and can be used in debug builds, where the
   compilation time is more important than the runtime.
   - Enable shared linking of ogs libraries. [#1133](https://github.com/ufz/ogs/pull/1133)
   - Break FileIO on ApplicationsLib dependency. [#1138](https://github.com/ufz/ogs/pull/1138)
   - Remove MeshLib on FileIO dependency. [#1143](https://github.com/ufz/ogs/pull/1143), [#1153](https://github.com/ufz/ogs/pull/1153)
   - Cleanup some of AssemblerLib dependencies. [#1166](https://github.com/ufz/ogs/pull/1166)
   - Split AssemblerLib and move to MathLib and NumLib [#1208](https://github.com/ufz/ogs/pull/1208)
   - Move InsituLib to MeshLib [#1208](https://github.com/ufz/ogs/pull/1208)
   - Remove MathLib depends on NumLib [#1208](https://github.com/ufz/ogs/pull/1208)
   - Remove dependency of FileIO on Data Explorer [#1302](https://github.com/ufz/ogs/pull/1302)
 - Introduced Conan package manager for automatic fetching of build dependencies, [#1141](https://github.com/ufz/ogs/pull/1141)
 - Inconsistent formatting of tabs and spaces was finally resolved: now all
   formatting, indentation and alignment, are done with four spaces. [#1201](https://github.com/ufz/ogs/pull/1201)
 - Windows 32-bit builds are disallowed because they are not supported.
   Can be forced by setting OGS_32_BIT=ON. [#1230](https://github.com/ufz/ogs/pull/1230)
 - Simplified FindEigen.cmake, [#1209](https://github.com/ufz/ogs/pull/1209)
 - git diff --check is run in its own Travis job, [#1207](https://github.com/ufz/ogs/pull/1207)

 - Moved some IO implementations from FileIO to BaseLib/IO, GeoLib/IO, MeshLib/IO, [#1182](https://github.com/ufz/ogs/pull/1182), [#1235](https://github.com/ufz/ogs/pull/1235)
 - Eigen is not optional anymore [#1218](https://github.com/ufz/ogs/pull/1218)
 - Removed OGS_USE_EIGENLIS CMake option. Use OGS_USE_LIS instead [#1251](https://github.com/ufz/ogs/pull/1251)

### Fixes:
 - Fix linking of Metis in MathLib. [#1147](https://github.com/ufz/ogs/pull/1147)
 - Fix memory leaks in GMSHInterface. [#1212](https://github.com/ufz/ogs/pull/1212)
 - Fix build with Lis [#1267](https://github.com/ufz/ogs/pull/1267)
 - Fixing several small issues with NetCDF import [#1169](https://github.com/ufz/ogs/pull/1169)
 - Restructure Applications related modules
    - Move DataHolderLib and FileIO under Applications [#1279](https://github.com/ufz/ogs/pull/1279)
 - Remove calling std::abort within libraries. Exeptions are thrown instead. [#1245](https://github.com/ufz/ogs/pull/1245)
 - Fix finding Boost with BOOST_ROOT CMake argument [#1287](https://github.com/ufz/ogs/pull/1287)
 - Fix linking of Sundials CVODE library [#1197](https://github.com/ufz/ogs/pull/1197)
 - Fixed issue where geometry names would be missing after merging geometries [#1295](https://github.com/ufz/ogs/pull/1295)

# 6.0.5

### Features:
 - Added an ODE solver library that can solve transient and nonlinear processes
   (see http://doxygen.opengeosys.org/df/d35/group__ODESolver.html).
 - Move up common Process parts from particular GroundwaterFlow process
   implementation. #951, #982
 - Separate Dirichlet boundary condition class implementation. #963
 - Split process output and post timestep. #998
 - Added pre- and postTimestep and -Iteration hooks to processes. [#1094](https://github.com/ufz/ogs/pull/1094), [#1100](https://github.com/ufz/ogs/pull/1100), [#1101](https://github.com/ufz/ogs/pull/1101)
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
   significantly. [#1092](https://github.com/ufz/ogs/pull/1092)
 - Command line argument `-l` for OGS cli and testrunner to specify verbosity of logging. [#1056](https://github.com/ufz/ogs/pull/1056)
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
- Optional support for VTK 7. [#1083](https://github.com/ufz/ogs/pull/1083)
- Test data is now a git submodule. [#1000](https://github.com/ufz/ogs/pull/1000)
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
