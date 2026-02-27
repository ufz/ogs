// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "ogs_embedded_python.h"

#include <pybind11/embed.h>

#include <filesystem>
#include <string>
#include <vector>

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

namespace
{
#ifndef _WIN32
std::vector<std::filesystem::path> findAlternativeSitePackagesPaths(
    std::filesystem::path const& venv_path)
{
    namespace fs = std::filesystem;

    std::vector<fs::path> alternatives;
    fs::path const lib_path = venv_path / "lib";

    if (!fs::exists(lib_path) || !fs::is_directory(lib_path))
    {
        // Should not happen, i.e. if venv directory is not valid
        return {};
    }

    for (auto const& entry : fs::directory_iterator(lib_path))
    {
        if (!entry.is_directory())
        {
            continue;
        }

        std::string const dirname = entry.path().filename().string();
        if (!dirname.starts_with("python"))
        {
            continue;
        }

        fs::path const candidate = entry.path() / "site-packages";
        if (fs::exists(candidate) && fs::is_directory(candidate))
        {
            alternatives.push_back(candidate);
        }
    }

    return alternatives;
}
#endif  // _WIN32

/// Finds site-packages path in the virtual environment
std::filesystem::path findSitePackagesPath(
    std::filesystem::path const& venv_path, int const emb_major,
    int const emb_minor)
{
    namespace fs = std::filesystem;

    // Construct path to site-packages directory, on *nix Python version is
    // embedded in the path. This is later used to check for compatibility.
#ifdef _WIN32
    // On Windows: venv/Lib/site-packages
    fs::path const site_packages = venv_path / "Lib" / "site-packages";
#else
    // On Unix: venv/lib/pythonX.Y/site-packages
    fs::path const site_packages = venv_path / "lib" /
                                   ("python" + std::to_string(emb_major) + "." +
                                    std::to_string(emb_minor)) /
                                   "site-packages";
#endif  // _WIN32

    if (!fs::exists(site_packages))
    {
#ifndef _WIN32
        // If correct site-packages directory is not found, check for
        // directories for other Python versions, indicating a Python version
        // mismatch.
        auto const alternatives = findAlternativeSitePackagesPaths(venv_path);
        if (!alternatives.empty())
        {
            std::string alternative_paths;
            for (std::size_t i = 0; i < alternatives.size(); ++i)
            {
                if (i > 0)
                {
                    alternative_paths += ", ";
                }
                alternative_paths += "'" + alternatives[i].string() + "'";
            }

            WARN(
                "Expected site-packages directory '{}' was not found. "
                "Found other site-packages directory/directories: {}. This "
                "may indicate a Python version mismatch between the embedded "
                "interpreter {}.{} and the virtual environment.",
                site_packages.string(), alternative_paths, emb_major,
                emb_minor);
        }
#endif  // _WIN32
        OGS_FATAL("site-packages directory not found at '{}'",
                  site_packages.string());
    }

    // Compatible venv directory found
    return site_packages;
}
}  // anonymous namespace

void setupEmbeddedPythonVenvPaths()
{
    namespace py = pybind11;
    namespace fs = std::filesystem;

    // Get embedded Python version
    py::object const version_info =
        py::module_::import("sys").attr("version_info");
    int const emb_major = version_info.attr("major").cast<int>();
    int const emb_minor = version_info.attr("minor").cast<int>();

    // Check for virtual environment
    char const* const venv = std::getenv("VIRTUAL_ENV");
    if (venv == nullptr)
    {
        DBUG("No virtual environment detected (VIRTUAL_ENV not set).");
        return;
    }

    fs::path const venv_path(venv);
    DBUG("Virtual environment detected at: {}", venv_path.string());

    // Find and validate site-packages path for the embedded interpreter
    // version. We do not execute the venv python binary because it might not
    // be runnable in all deployment setups.
    fs::path const site_packages =
        findSitePackagesPath(venv_path, emb_major, emb_minor);
    INFO("Using virtual environment site-packages: {}", site_packages.string());

    // Add to sys.path (insert at beginning for highest priority)
    py::list sys_path = py::module_::import("sys").attr("path");
    sys_path.insert(0, py::str(site_packages.string()));
}

}  // namespace ApplicationsLib
