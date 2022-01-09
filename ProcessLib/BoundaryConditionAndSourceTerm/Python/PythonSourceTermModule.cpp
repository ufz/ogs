/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PythonSourceTermModule.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "PythonSourceTermPythonSideInterface.h"

namespace ProcessLib
{
namespace SourceTerms
{
namespace Python
{
//! Trampoline class allowing methods of class
//! PythonSourceTermPythonSideInterface to be overridden on the Python
//! side. Cf. https://pybind11.readthedocs.io/en/stable/advanced/classes.html
class PythonSourceTermPythonSideInterfaceTrampoline
    : public PythonSourceTermPythonSideInterface
{
public:
    using PythonSourceTermPythonSideInterface::
        PythonSourceTermPythonSideInterface;

    std::pair<double, std::vector<double>> getFlux(
        double t, std::array<double, 3> const& x,
        std::vector<double> const& primary_variables) const override
    {
        using Ret = std::pair<double, std::vector<double>>;
        PYBIND11_OVERLOAD_PURE(Ret, PythonSourceTermPythonSideInterface,
                               getFlux, t, x, primary_variables);
    }
};

void pythonBindSourceTerm(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<PythonSourceTermPythonSideInterface,
               PythonSourceTermPythonSideInterfaceTrampoline>
        pybc(m, "SourceTerm");

    pybc.def(py::init());

    pybc.def("getFlux", &PythonSourceTermPythonSideInterface::getFlux);
}

}  // namespace Python
}  // namespace SourceTerms
}  // namespace ProcessLib
