/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifdef OGS_USE_PYTHON

#include <pybind11/embed.h>
#include <pybind11/stl.h>
#include <logog/include/logog.hpp>

#include "PythonBoundaryConditionDetail.h"

namespace ProcessLib
{
void pythonBindBoundaryCondition(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<PythonBoundaryConditionPythonSideInterface,
               PythonBoundaryConditionPythonSideInterfaceTrampoline>
        pybc(m, "BoundaryCondition");

    pybc.def(py::init());

    pybc.def("getDirichletBCValue",
             &PythonBoundaryConditionPythonSideInterface::getDirichletBCValue);
    pybc.def("getFlux", &PythonBoundaryConditionPythonSideInterface::getFlux);
}

}  // namespace ProcessLib

PYBIND11_EMBEDDED_MODULE(OpenGeoSys, m)
{
    DBUG("Binding Python module OpenGeoSys.");

    ProcessLib::pythonBindBoundaryCondition(m);
}

#endif  // OGS_USE_PYTHON
