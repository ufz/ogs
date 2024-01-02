/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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

    std::tuple<double, std::vector<double>, std::vector<double>,
               std::vector<int>, std::vector<double>>
    initializeDataContainer() const override
    {
        using Ret = std::tuple<double, std::vector<double>, std::vector<double>,
                               std::vector<int>, std::vector<double>>;
        PYBIND11_OVERLOAD(Ret,
                          BHEInflowPythonBoundaryConditionPythonSideInterface,
                          initializeDataContainer);
    }

    std::tuple<bool, bool, std::vector<double>, std::vector<double>>
    tespySolver(double t,
                std::vector<double> const& Tin_val,
                std::vector<double> const& Tout_val) const override
    {
        using Ret =
            std::tuple<bool, bool, std::vector<double>, std::vector<double>>;
        PYBIND11_OVERLOAD(Ret,
                          BHEInflowPythonBoundaryConditionPythonSideInterface,
                          tespySolver, t, Tin_val, Tout_val);
    }

    std::tuple<std::vector<double>, std::vector<double>>
    serverCommunicationPreTimestep(
        double const t, double const dt, std::vector<double> const& Tin_val,
        std::vector<double> const& Tout_val,
        std::vector<double> const& flowrate) const override
    {
        using Ret = std::tuple<std::vector<double>, std::vector<double>>;
        PYBIND11_OVERLOAD(
            Ret, BHEInflowPythonBoundaryConditionPythonSideInterface,
            serverCommunicationPreTimestep, t, dt, Tin_val, Tout_val, flowrate);
    }

    void serverCommunicationPostTimestep(
        double const t, double const dt, std::vector<double> const& Tin_val,
        std::vector<double> const& Tout_val,
        std::vector<double> const& flowrate) const override
    {
        PYBIND11_OVERLOAD(void,
                          BHEInflowPythonBoundaryConditionPythonSideInterface,
                          serverCommunicationPostTimestep, t, dt, Tin_val,
                          Tout_val, flowrate);
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
             &BHEInflowPythonBoundaryConditionPythonSideInterface::
                 initializeDataContainer);
    pybc.def("tespySolver",
             &BHEInflowPythonBoundaryConditionPythonSideInterface::tespySolver);

    pybc.def("serverCommunicationPreTimestep",
             &BHEInflowPythonBoundaryConditionPythonSideInterface::
                 serverCommunicationPreTimestep);

    pybc.def("serverCommunicationPostTimestep",
             &BHEInflowPythonBoundaryConditionPythonSideInterface::
                 serverCommunicationPostTimestep);
}

}  // namespace ProcessLib
