# 6.0.2

| Release planned 2015/06/15, [GitHub Release Link](https://github.com/ufz/ogs/releases/tag/6.0.2)

## Release notes

The second release of ogs6 introduces Neumann boundary conditions and VTK result output.

### Features:

- [Neumann boundary conditions](http://docs.opengeosys.org/docs/benchmarks/elliptic/groundwater-flow-neumann)
- Implement mesh properties for storage of data fields. This also includes mapped values (e.g. based on material id) [PR #542](https://github.com/ufz/ogs/pull/542), [PR #644](https://github.com/ufz/ogs/pull/644), 
- Refactored mesh property classes to [enable VTK output](https://github.com/ufz/ogs/pull/692)
- Extended the available elements to quadratic (e.g. Quad8) based on generalized "element rules".  [#572](https://github.com/ufz/ogs/pull/572), [PR #656](https://github.com/ufz/ogs/pull/656), [PR #657](https://github.com/ufz/ogs/pull/657),
- Extend computation of shape matrices to lower dimensional elements embedded in higher dimensional space [#655]((https://github.com/ufz/ogs/pull/655)
- Builds with MinGW (GCC) on Windows, see [Developer Guide](http://docs.opengeosys.org/docs/devguide/getting-started/prerequisites) and the new MinGW platform instructions
- [Cross-compiling for Windows with MXE on Mac OS](http://docs.opengeosys.org/docs/devguide/advanced/cross-compiling)
- [Support for new cross-platform IDE CLion](http://docs.opengeosys.org/docs/devguide/development-workflows/development-ides#clion)
- Add gradual refinement to the structured mesh generator [PR #539](https://github.com/ufz/ogs/pull/539)
- Add a command line tool "queryMesh" to search mesh information [PR #665](https://github.com/ufz/ogs/pull/665)

### Fixes

- Performance optimizations in VTK mesh conversion, [PR #695](https://github.com/ufz/ogs/pull/695)
- Improve layered prism mesh construction
- Fix a lot of warnings on gcc/clang/msvc compilers improving the code


### Infrastructure

- Test runtime is monitored at Jenkins for [normal](https://svn.ufz.de:8443/job/OGS-6/job/Tests_Linux/) and [nightly large tests](https://svn.ufz.de:8443/job/OGS-6/job/Tests_Linux_Large/)
- Utilities are build by separate Jenkins Jobs, e.g. [Win_Utils](https://svn.ufz.de:8443/job/OGS-6/job/Win_Utils/)

## Test examples
![](http://docs.opengeosys.org/assets/files/square_1e2_neumann_result.png)

 - Example 1: TODO

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
