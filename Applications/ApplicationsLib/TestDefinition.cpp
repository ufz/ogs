/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TestDefinition.h"

#include <cstdlib>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "BaseLib/FileTools.h"
#ifdef USE_PETSC
#include "MeshLib/IO/VtkIO/VtuInterface.h"  // For petsc file name conversion.
#include <petsc.h>
#endif

namespace
{
/// Safe conversion of a double to a string using decimal or decimal exponent
/// notation. See std::snprintf() for details for "%g" conversion specifier.
/// \note std::to_string uses "%f" conversion specifier which is not sufficient
/// for small numbers like 1e-15.
std::string convert_to_string(double const& value)
{
    // TODO (naumov) Replace this with fmt library.
    char buffer[32];
    int const chars_written =
        std::snprintf(buffer, sizeof(buffer), "%g", value);
    if (chars_written < 0 || chars_written >= static_cast<int>(sizeof(buffer)))
    {
        OGS_FATAL("Could not convert a value to string.");
    }
    return std::string{buffer};
}

/// Wraps a string into double ticks.
std::string safeString(std::string s)
{
    return "\"" + s + "\"";
}

/// Tries to find a vtkdiff executable by testing 'path/vtkdiff --version' calls
/// for various paths.
std::string findVtkdiff()
{
    // Try to read the VTKDIFF_EXE environment variable.
    if (const char* vtkdiff_exe_environment_variable =
            std::getenv("VTKDIFF_EXE"))
    {
        std::string const vtkdiff_exe{vtkdiff_exe_environment_variable};
        DBUG("VTKDIFF_EXE set to %s.", vtkdiff_exe.c_str());

        //
        // Sanity checks.
        //
        {  // The base name is 'vtkdiff'
            auto const& base_name =
                BaseLib::extractBaseNameWithoutExtension(vtkdiff_exe);
            if (base_name != "vtkdiff")
            {
                OGS_FATAL(
                    "The VTKDIFF_EXE environment variable does not point to "
                    "'vtkdiff'. VTKDIFF_EXE='%s'",
                    vtkdiff_exe.c_str());
            }
        }
        {  // vtkdiff must exist.
            if (!BaseLib::IsFileExisting(vtkdiff_exe))
            {
                OGS_FATAL(
                    "The VTKDIFF_EXE points to a non-existing file. "
                    "VTKDIFF_EXE='%s'",
                    vtkdiff_exe.c_str());
            }
        }

        //
        // Test the actual call.
        //
        int const return_value =
            // TODO (naumov) replace system call with output consuming call
            // (fork + execl seems to be more safe), and extract the vtkdiff
            // call to common function. Also properly escape all strings in
            // command lines.
            // Reference for POSIX and Windows:
            // https://wiki.sei.cmu.edu/confluence/pages/viewpage.action?pageId=87152177
            // Take care when using fork, which might copy resources.
            std::system((vtkdiff_exe + " --version").c_str());
        if (return_value == 0)
        {
            return vtkdiff_exe;
        }
        WARN(
            "Calling %s from the VTKDIFF_EXE environment variable didn't work "
            "as expected. Return value was %d.",
            vtkdiff_exe.c_str(), return_value);
    }

    std::string const vtkdiff_exe{"vtkdiff"};
    std::vector<std::string> const paths = {"", "bin"};
    auto const path = find_if(
        begin(paths), end(paths), [&vtkdiff_exe](std::string const& path) {
            int const return_value =
                // TODO (naumov) replace system call with output consuming call
                // as in an above todo comment.
                std::system(
                    (BaseLib::joinPaths(path, vtkdiff_exe) + " --version")
                        .c_str());
            return return_value == 0;
        });
    if (path == end(paths))
    {
        OGS_FATAL("vtkdiff not found.");
    }
    return BaseLib::joinPaths(*path, vtkdiff_exe);
}

}  // namespace

namespace ApplicationsLib
{
TestDefinition::TestDefinition(BaseLib::ConfigTree const& config_tree,
                               std::string const& reference_path,
                               std::string const& output_directory)
{
    if (reference_path.empty())
    {
        OGS_FATAL(
            "Reference path containing expected result files can not be "
            "empty.");
    }

    std::string const vtkdiff = findVtkdiff();

    // Construct command lines for each entry.
    auto const& vtkdiff_configs = config_tree.getConfigSubtreeList("vtkdiff");
    _command_lines.reserve(vtkdiff_configs.size());
    for (auto const& vtkdiff_config : vtkdiff_configs)
    {
        std::string const& field_name =
            vtkdiff_config.getConfigParameter<std::string>("field");
        DBUG("vtkdiff will compare field '%s'.", field_name.c_str());

#ifdef USE_PETSC
        int rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        int mpi_size;
        MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
        std::string const& filename =
            MeshLib::IO::getVtuFileNameForPetscOutputWithoutExtension(
                vtkdiff_config.getConfigParameter<std::string>("file")) +
            "_" + std::to_string(rank) + ".vtu";
#else
        std::string const& filename =
            vtkdiff_config.getConfigParameter<std::string>("file");
#endif  // OGS_USE_PETSC
        std::string const& output_filename =
            BaseLib::joinPaths(output_directory, filename);
        _output_files.push_back(output_filename);
        // TODO (naumov) expand filename relative to ref path for globbing.
        std::string const& reference_filename =
            BaseLib::joinPaths(reference_path, filename);

        auto const& absolute_tolerance =
            vtkdiff_config.getConfigParameterOptional<double>(
                "absolute_tolerance");
        std::string const absolute_tolerance_parameter =
            absolute_tolerance == boost::none
                ? ""
                : "--abs " + convert_to_string(*absolute_tolerance);

        auto const& relative_tolerance =
            vtkdiff_config.getConfigParameterOptional<double>(
                "relative_tolerance");
        std::string const relative_tolerance_parameter =
            relative_tolerance == boost::none
                ? ""
                : "--rel " + convert_to_string(*relative_tolerance);

        //
        // Construct command line.
        //
        std::string command_line =
            vtkdiff + " -a " + safeString(field_name) + " -b " +
            safeString(field_name) + " " + safeString(reference_filename) +
            " " + safeString(output_filename) + " " +
            absolute_tolerance_parameter + " " + relative_tolerance_parameter;
        INFO("Will run '%s'", command_line.c_str());
        _command_lines.emplace_back(std::move(command_line));
    }
}

bool TestDefinition::runTests() const
{
    std::vector<int> return_values;
    transform(begin(_command_lines), end(_command_lines),
              back_inserter(return_values),
              [](std::string const& command_line) {
                  int const return_value = std::system(command_line.c_str());
                  if (return_value != 0)
                  {
                      WARN("Return value %d was returned by '%s'.",
                           return_value, command_line.c_str());
                  }
                  return return_value;
              });
    return !return_values.empty() &&
           all_of(begin(return_values), end(return_values),
                  [](int const& return_value) { return return_value == 0; });
}

std::vector<std::string> const& TestDefinition::getOutputFiles() const
{
    return _output_files;
}
}  // namespace ApplicationsLib
