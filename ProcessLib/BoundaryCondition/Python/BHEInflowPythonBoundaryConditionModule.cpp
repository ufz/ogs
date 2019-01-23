/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BHEInflowPythonBoundaryConditionModule.h"

#include <pybind11/stl.h>

#include "BHEInflowPythonBoundaryConditionPythonSideInterface.h"

namespace ProcessLib
{
//! Trampoline class allowing methods of class
//! PythonBoundaryConditionPythonSideInterface to be overridden on the Python
//! side. Cf. https://pybind11.readthedocs.io/en/stable/advanced/classes.html
class BHEInflowPythonBoundaryConditionPythonSideInterfaceTrampoline
    : public BHEInflowPythonBoundaryConditionPythonSideInterface
{
public:
    using BHEInflowPythonBoundaryConditionPythonSideInterface::
        BHEInflowPythonBoundaryConditionPythonSideInterface;

    std::tuple<bool, std::vector<double>, std::vector<double>, std::vector<int>, double>
    initializeDataContainer() const override
    {
        using Ret = std::tuple<bool, std::vector<double>, std::vector<double>, std::vector<int>, double>;
        PYBIND11_OVERLOAD(Ret,
                          BHEInflowPythonBoundaryConditionPythonSideInterface,
                          initializeDataContainer);
    }

    std::tuple<bool, bool, std::vector<double>> tespyThermalSolver(
        double t,
        std::vector<double> const& Tin_val,
        std::vector<double> const& Tout_val) const override
    {
        using Ret = std::tuple<bool, bool, std::vector<double>>;
        PYBIND11_OVERLOAD(Ret, BHEInflowPythonBoundaryConditionPythonSideInterface,
                          tespyThermalSolver, t, Tin_val, Tout_val);
    }

    std::tuple<bool, std::vector<double>> tespyHydroSolver() const override
    {
        using Ret = std::tuple<bool, std::vector<double>>;
        PYBIND11_OVERLOAD(Ret,
                          BHEInflowPythonBoundaryConditionPythonSideInterface,
                          tespyHydroSolver);
    }
};

void bheInflowpythonBindBoundaryCondition(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<BHEInflowPythonBoundaryConditionPythonSideInterface,
               BHEInflowPythonBoundaryConditionPythonSideInterfaceTrampoline>
        pybc(m, "BHENetwork");

    pybc.def(py::init());

    pybc.def("initializeDataContainer",
        &BHEInflowPythonBoundaryConditionPythonSideInterface::initializeDataContainer);
    pybc.def("tespyThermalSolver",
        &BHEInflowPythonBoundaryConditionPythonSideInterface::tespyThermalSolver);
    pybc.def("tespyHydroSolver",
             &BHEInflowPythonBoundaryConditionPythonSideInterface::tespyHydroSolver);
}

}  // namespace ProcessLib
