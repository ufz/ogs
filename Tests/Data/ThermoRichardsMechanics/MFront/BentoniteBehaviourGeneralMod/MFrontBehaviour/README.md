# Bentonite Material Behaviour implemented in generalmod used in OGS via MFront

This directory contains the implementation of a bentonite material behaviour
implemented in generalmod from the TRIAX Element Test Driver by
David Mašín et al.
TRIAX can be downloaded from https://soilmodels.com/download/triax-zip/.

Next to the generalmod files this directory contains the associated MFront
behaviour definition and a header file converting data between MFront and
generalmod.

The material parameters in the MFront file have been calibrated by D. Mašín to a
number of experiments. Therefore, these parameters (hence, the entire MFront
behaviour) are material specific and we don't include it in OGS's library of
shipped MFront models, but use it for some tests only.

## Copyright/License

generalmod is licensed under a 3 clause BSD-license, Copyright (c) David Mašín <masin@natur.cuni.cz>.
The other files in this directory are subject to OpenGeoSys's usual 3-clause BSD license.

## Acknowledgements

The following people were involved in the implementation and testing that
finally led to successfully interfacing generalmod with OpenGeoSys via MFront:

Éric Simo, Thomas Nagel, Thomas Helfer, David Mašín, Tymofiy Gerasimov,
Tomáš Krejčí, Christoph Lehmann.

## References

* Mašín, David. 2013. “Clay Hypoplasticity with Explicitly Defined Asymptotic
  States.” Acta Geotechnica 8: 481–96.
* Mašín, David. 2017. “Coupled Thermohydromechanical Double-Structure Model for
  Expansive Soils.” Journal of Engineering Mechanics 143 (9): 04017067.
  https://doi.org/10.1061/(ASCE)EM.1943-7889.0001278.
