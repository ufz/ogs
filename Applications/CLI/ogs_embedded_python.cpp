/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ogs_embedded_python.h"

#include <pybind11/embed.h>

#include <algorithm>

#include "BaseLib/Logging.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/Python/BHEInflowPythonBoundaryConditionModule.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/Python/PythonBoundaryConditionModule.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/Python/PythonSourceTermModule.h"

PYBIND11_EMBEDDED_MODULE(OpenGeoSys, m)
{
    DBUG("Binding Python module OpenGeoSys.");

    ProcessLib::pythonBindBoundaryCondition(m);
    ProcessLib::bheInflowpythonBindBoundaryCondition(m);
    ProcessLib::SourceTerms::Python::pythonBindSourceTerm(m);

    // Check for activated virtual environment and add it to sys.path
    pybind11::exec(R"(
        import os
        import sys
        if "VIRTUAL_ENV" in os.environ:
            venv_site_packages_path = f"{os.environ['VIRTUAL_ENV']}/lib/python{sys.version_info.major}.{sys.version_info.minor}/site-packages"
            if os.path.exists(venv_site_packages_path):
                print(
                    f"Virtual environment detected, adding {venv_site_packages_path} to sys.path."
                )
                sys.path.insert(0, venv_site_packages_path)
    )");
}

namespace ApplicationsLib
{
pybind11::scoped_interpreter setupEmbeddedPython()
{
    // Allows ogs to be interrupted by SIGINT, which otherwise is handled by
    // python. See
    // https://docs.python.org/3/c-api/exceptions.html#c.PyErr_CheckSignals and
    // https://pybind11.readthedocs.io/en/stable/faq.html#how-can-i-properly-handle-ctrl-c-in-long-running-functions
    constexpr bool init_signal_handlers = false;
    return pybind11::scoped_interpreter{init_signal_handlers};
}

}  // namespace ApplicationsLib
