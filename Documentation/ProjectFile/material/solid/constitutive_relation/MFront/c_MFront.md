Loads a material model from a shared library that has been created by MFront
with MFront's generic interface.

See <https://thelfer.github.io/tfel/web/> <https://github.com/thelfer/MFrontGenericInterfaceSupport>

This requires MFront version 3.2.1 (which must correspond to the MFront's
"generic" interface), which must be available on your system: see
<https://thelfer.github.io/tfel/web> page for download and installation instructions.

\note The MFront library is distributed under GPL license.

\attention
OpenGeoSys initializes the internal state of the MFront material model to
all-zero. Furthermore, the very first iteration of OpenGeoSys is performed
already at the first load step, i.e., OpenGeoSys does not perform a stress
integration with the IC values. That differs, e.g., from the `mtest` executable
shipped with MFront, which, initially does a stress integration with the ICs
values. If your model needs such a stress integration step, you must make sure
that you prepend an initial, artificial load step. If you observe convergence
issues during the local stress integration in the very first iteration of
OpenGeoSys, you should consider adding such a load step, in particular if `mtest`
succeeds to integrate the same load path.
