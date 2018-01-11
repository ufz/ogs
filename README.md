OpenGeoSys 6
============

[![Tag](https://img.shields.io/github/tag/ufz/ogs.svg?style=flat-square)](https://github.com/ufz/ogs/releases)
[![BSD License (modified)](http://img.shields.io/badge/license-BSD-blue.svg?style=flat-square)](https://github.com/ufz/ogs/blob/master/LICENSE.txt)
[![Build Status](https://jenkins.opengeosys.org/buildStatus/icon?job=ufz/ogs/master)](https://jenkins.opengeosys.org/job/ufz/job/ogs/job/master)
[![DOI](https://zenodo.org/badge/1701384.svg)](https://zenodo.org/badge/latestdoi/1701384)

[OpenGeoSys][ogs] (OGS) is a scientific open source project for the development of
numerical methods for the simulation of thermo-hydro-mechanical-chemical
(THMC) processes in porous and fractured media. OGS is implemented in C++, it
is object-oriented with an focus on the numerical solution of coupled multi-field
problems (multi-physics). Parallel versions of OGS are available relying on
both MPI and OpenMP concepts. Application areas of OGS are currently CO2
sequestration, geothermal energy, water resources management, hydrology and
waste deposition. OGS is comprised of the THMC-simulator (simply referred to as
*OGS*) and a visualization tool (*Data Explorer*). OGS is developed by the
[OpenGeoSys Community][ogs].

## Current status ##

It has to be noticed that the current OGS-6 version is under very heavy development
and [vivid discussion](https://github.com/ufz/ogs/issues), and does not implement all
functionality from OGS-5.

## Software development ##

- Good starting point for users as well as for developers is the [documentation][documentation]
- Check your code against our [styleguide](http://ufz.github.io/styleguide/cppguide.xml)
- Have a look at the [source code documentation][docs]
- For the actual build status see the [Jenkins-CI server][jenkins-ci]

## License ##

OpenGeoSys is distributed under a Modified BSD License which encourages users to
attribute the work of the OpenGeoSys Community especially in scientific
publications. See the [LICENSE.txt][license-source] for the license text.

[ogs]: http://www.opengeosys.org
[documentation]: https://docs.opengeosys.org/docs/
[jenkins-ci]: https://jenkins.opengeosys.org/job/ufz/job/ogs/job/master/
[docs]: http://doxygen.opengeosys.org
[license-source]: https://github.com/ufz/ogs/blob/master/LICENSE.txt
