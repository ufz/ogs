// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "ogs_embedded_python.h"

#include <pybind11/embed.h>

#include <cstdio>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>

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
/// Custom deleter for FILE handles from popen
struct PipeCloser
{
#ifdef _WIN32
    void operator()(FILE* f) const { _pclose(f); }
#else
    void operator()(FILE* f) const { pclose(f); }
#endif  // _WIN32
};

/// Executes a command and captures its stdout output using popen
std::optional<std::string> executeCommand(std::string_view command)
{
    std::array<char, 256> buffer;
    std::string result;

#ifdef _WIN32
    std::unique_ptr<FILE, PipeCloser> pipe(_popen(command.data(), "r"));
#else
    std::unique_ptr<FILE, PipeCloser> pipe(popen(command.data(), "r"));
#endif  // _WIN32

    if (!pipe)
    {
        DBUG("Failed to execute command: {}", command);
        return std::nullopt;
    }

    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
    {
        result += buffer.data();
    }

    return result;
}

/// Gets Python version by executing the venv's python executable
std::optional<std::pair<int, int>> getPythonVersionFromVenv(
    std::filesystem::path const& venv_path)
{
    namespace fs = std::filesystem;

    fs::path python_exe;
#ifdef _WIN32
    python_exe = venv_path / "Scripts" / "python.exe";
#else
    python_exe = venv_path / "bin" / "python";
#endif  // _WIN32

    if (!fs::exists(python_exe))
    {
        DBUG("Python executable not found at: {}", python_exe.string());
        return std::nullopt;
    }

    std::string const command = python_exe.string() + " --version";
    auto const output = executeCommand(command);

    if (!output.has_value())
    {
        DBUG("Failed to get Python version from: {}", python_exe.string());
        return std::nullopt;
    }

    // Parse output like "Python 3.11.5"
    std::string_view const out_view(output.value());
    constexpr std::string_view prefix = "Python ";
    if (!out_view.starts_with(prefix))
    {
        DBUG("Unexpected Python version output: {}", output.value());
        return std::nullopt;
    }

    std::string_view version_part = out_view.substr(prefix.size());
    int major = 0, minor = 0;

    if (std::sscanf(version_part.data(), "%d.%d", &major, &minor) != 2)
    {
        DBUG("Failed to parse Python version from: {}", output.value());
        return std::nullopt;
    }

    DBUG("Detected Python version {}.{} from venv", major, minor);
    return std::make_pair(major, minor);
}

/// Finds site-packages path in the virtual environment
std::filesystem::path findSitePackagesPath(
    std::filesystem::path const& venv_path)
{
    namespace fs = std::filesystem;

#ifdef _WIN32
    // On Windows: venv/Lib/site-packages
    fs::path const site_packages = venv_path / "Lib" / "site-packages";
#else
    // On Unix: venv/lib/pythonX.Y/site-packages
    auto const version = getPythonVersionFromVenv(venv_path);
    if (!version.has_value())
    {
        OGS_FATAL(
            "Failed to determine Python version of the virtual environment.");
    }
    fs::path const site_packages = venv_path / "lib" /
                                   ("python" + std::to_string(version->first) +
                                    "." + std::to_string(version->second)) /
                                   "site-packages";
#endif  // _WIN32

    if (!fs::exists(site_packages))
    {
        OGS_FATAL("site-packages directory not found at '{}'",
                  site_packages.string());
    }

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

    auto const venv_version = getPythonVersionFromVenv(venv_path);
    if (!venv_version.has_value())
    {
        OGS_FATAL(
            "Failed to determine Python version from virtual environment at "
            "'{}'. "
            "Please ensure the virtual environment is valid.",
            venv_path.string());
    }

    int const venv_major = venv_version->first;
    int const venv_minor = venv_version->second;

    // Validate version match
    if (venv_major != emb_major || venv_minor != emb_minor)
    {
        OGS_FATAL(
            "Python version mismatch:\n"
            "  Embedded interpreter: {}.{}\n"
            "  Virtual environment:  {}.{}\n"
            "The virtual environment must use the same Python version as the "
            "embedded interpreter.",
            emb_major, emb_minor, venv_major, venv_minor);
    }

    // Find and validate site-packages path
    fs::path const site_packages = findSitePackagesPath(venv_path);
    INFO("Using virtual environment site-packages: {}", site_packages.string());

    // Add to sys.path (insert at beginning for highest priority)
    py::list sys_path = py::module_::import("sys").attr("path");
    sys_path.insert(0, py::str(site_packages.string()));
}

}  // namespace ApplicationsLib
