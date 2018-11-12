Loads a material model from a shared library that has been created by MFront
with MFront's generic interface.

See http://tfel.sourceforge.net/ https://github.com/thelfer/MFrontGenericInterfaceSupport

This requires MFront 3.2.0 or newer (because of MFront's "generic" interface).

\attention
OpenGeoSys initializes the internal state of the MFront material model to
all-zero. Furthermore, the very first iteration of OpenGeoSys is performed
already at the first load step, i.e., OpenGeoSys does not perform a stress
integration with the IC values. That differs, e.g., from the mtest executable
shipped with MFront, which, initially does a stress integration with the ICs
values. If your model needs such a stress integration step, you must make sure
that you prepend an initial, artificial load step. If you observe convergence
issues during the local stress integration in the very first iteration of
OpenGeoSys, you should consider adding such a load step, in particular if mtest
succeeds to integrate the same load path.
