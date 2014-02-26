OpenGeoSys 6                                                      {#mainpage}
============

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

It has to be noticed that the current ogs6 version is under very heavy development
and [vivid discussion](https://github.com/ufz/ogs/issues), and not even reaching a
beta release version, i.e. not all data structure and features are fixed yet. If
you would like to have a peek into the current status of the code, we recommend to
have a look at Norihiro's [draft branch](https://github.com/norihiro-w/ogs/tree/draft).
Which is merely a proposed prototype. However, it should compile out of the box
and supports many of the simulation processes in ogs5 already. The plan is to
refactor the code and to merge in functionality into the official ogs6-repository
feature-by-feature.

## Software development ##

- To get started checkout the [developer guide][devguide]
- Check your code against our [styleguide](http://ufz.github.io/styleguide/cppguide.xml)
- Have a look at the [source code documentation][docs]
- For the actual build status see the [Jenkins-CI server][jenkins-ci]
- Actual Travis build status: [![Build Status](https://travis-ci.org/ufz/ogs.png)](https://travis-ci.org/ufz/ogs)

## License ##

OpenGeoSys is distributed under a Modified BSD License which encourages users to
attribute the work of the OpenGeoSys Community especially in scientific
publications. See the [LICENSE.txt][license-source] for the license text.

[ogs]: http://www.opengeosys.com
[devguide]: http://docs.opengeosys.org/docs/devguide
[jenkins-ci]: https://svn.ufz.de/jenkins/job/OGS-6/
[docs]: https://svn.ufz.de/jenkins/job/OGS-6/job/Docs/lastSuccessfulBuild/artifact/build/docs/index.html
[license-source]: https://github.com/ufz/ogs/blob/master/LICENSE.txt
