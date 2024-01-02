/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <pybind11/pybind11.h>

#include "ProcessLib/BoundaryConditionAndSourceTerm/Python/BHEInflowPythonBoundaryConditionModule.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/Python/PythonBoundaryConditionModule.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/Python/PythonSourceTermModule.h"

PYBIND11_MODULE(callbacks, m)
{
    m.attr("__name__") = "ogs.callbacks";
    ProcessLib::pythonBindBoundaryCondition(m);
    ProcessLib::bheInflowpythonBindBoundaryCondition(m);
    ProcessLib::SourceTerms::Python::pythonBindSourceTerm(m);
}
