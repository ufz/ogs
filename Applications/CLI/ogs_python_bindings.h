#pragma once

#ifdef OGS_USE_PYTHON

#include <pybind11/embed.h>

namespace ApplicationsLib
{
/// Sets up an embedded Python interpreter and makes sure that the OpenGeoSys
/// Python module is not removed by the linker.
pybind11::scoped_interpreter setupEmbeddedPython();

}  // namespace ApplicationsLib

#endif  // OGS_USE_PYTHON
