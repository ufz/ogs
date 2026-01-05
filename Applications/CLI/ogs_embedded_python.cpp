// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "ogs_embedded_python.h"

#include <pybind11/embed.h>

#include <algorithm>
#include <filesystem>

#include "BaseLib/Error.h"
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

void setupEmbeddedPythonVenvPaths()
{
    namespace py = pybind11;
    namespace fs = std::filesystem;

    py::module_ sys = py::module_::import("sys");
    py::list sys_path = sys.attr("path");

    int emb_major = sys.attr("version_info").attr("major").cast<int>();
    int emb_minor = sys.attr("version_info").attr("minor").cast<int>();

    char const* venv = std::getenv("VIRTUAL_ENV");
    if (!venv)
    {
        DBUG("NO virtual environment detected.");
        return;
    }

    fs::path venv_lib = fs::path(venv) / "lib";
    if (!fs::exists(venv_lib))
    {
        DBUG("VIRTUAL_ENV set but {} does not exist.", venv_lib.string());
        return;
    }

    std::optional<std::pair<int, int>> venv_version;
    fs::path python_dir;

    for (auto const& entry : fs::directory_iterator(venv_lib))
    {
        if (!entry.is_directory())
            continue;

        auto name = entry.path().filename().string();
        int maj = 0, min = 0;

        if (std::sscanf(name.c_str(), "python%d.%d", &maj, &min) == 2)
        {
            venv_version = {maj, min};
            python_dir = entry.path();
            break;
        }
    }

    if (!venv_version)
    {
        OGS_FATAL(
            "Virtual environment at {} does not contain a pythonX.Y directory.",
            venv_lib.string());
    }

    auto [venv_major, venv_minor] = *venv_version;

    if (venv_major != emb_major || venv_minor != emb_minor)
    {
        OGS_FATAL(
            "Python version mismatch:\n"
            "  Embedded interpreter : {}.{}\n"
            "  Virtual environment  : {}.{}\n"
            "The virtual environment must be created with the same Python "
            "version as OGS was compiled with.",
            emb_major, emb_minor, venv_major, venv_minor);
    }

    fs::path site_packages = python_dir / "site-packages";
    if (!fs::exists(site_packages))
    {
        OGS_FATAL("site-packages directory not found at {}",
                  site_packages.string());
    }

    DBUG("Using virtual environment site-packages: {}", site_packages.string());

    sys_path.insert(0, py::str(site_packages.string()));
}

}  // namespace ApplicationsLib
