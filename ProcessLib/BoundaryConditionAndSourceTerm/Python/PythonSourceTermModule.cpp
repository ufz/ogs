// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "PythonSourceTermModule.h"

#include <algorithm>
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
