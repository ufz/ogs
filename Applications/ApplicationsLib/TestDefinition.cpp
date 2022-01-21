/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TestDefinition.h"

#include <cmath>
#include <cstdlib>
#include <regex>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "BaseLib/FileTools.h"
#include "filesystem.h"
#ifdef USE_PETSC
#include <petsc.h>

#include "MeshLib/IO/VtkIO/VtuInterface.h"  // For petsc file name conversion.
#endif

namespace
{
/// Test if the given string is convertible to a valid double value, not a NaN.
bool isConvertibleToDouble(std::string const& s)
{
    std::size_t pos = 0;
    double value;
    try
    {
        value = std::stod(s, &pos);
    }
    catch (...)
    {
        OGS_FATAL("The given string '{:s}' is not convertible to double.", s);
    }
    if (pos != s.size())
    {
        OGS_FATAL(
            "Only {:d} characters were used for double conversion of string "
            "'{:s}'",
            pos, s);
    }

    if (std::isnan(value))
    {
        OGS_FATAL("The given string '{:s}' results in a NaN value.", s);
    }
    return true;
}

/// Wraps a string into double ticks.
std::string safeString(std::string const& s)
{
    std::stringstream ss;
    ss << std::quoted(s);
    return ss.str();
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
        DBUG("VTKDIFF_EXE set to {:s}.", vtkdiff_exe);

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
                    "'vtkdiff'. VTKDIFF_EXE='{:s}'",
                    vtkdiff_exe);
            }
        }
        {  // vtkdiff must exist.
            if (!BaseLib::IsFileExisting(vtkdiff_exe))
            {
                OGS_FATAL(
                    "The VTKDIFF_EXE points to a non-existing file. "
                    "VTKDIFF_EXE='{:s}'",
                    vtkdiff_exe);
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
            "Calling {:s} from the VTKDIFF_EXE environment variable didn't "
            "work as expected. Return value was {:d}.",
            vtkdiff_exe, return_value);
    }

    std::string const vtkdiff_exe{"vtkdiff"};
    std::vector<std::string> const paths = {"", "bin"};
    auto const path =
        find_if(begin(paths), end(paths),
                [&vtkdiff_exe](std::string const& path)
                {
                    int const return_value =
                        // TODO (naumov) replace system call with output
                        // consuming call as in an above todo comment.
                        std::system((BaseLib::joinPaths(path, vtkdiff_exe) +
                                     " --version")
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
    //! \ogs_file_param{prj__test_definition__vtkdiff}
    auto const& vtkdiff_configs = config_tree.getConfigSubtreeList("vtkdiff");
    _command_lines.reserve(vtkdiff_configs.size());
    for (auto const& vtkdiff_config : vtkdiff_configs)
    {
        std::string const& field_name =
            //! \ogs_file_param{prj__test_definition__vtkdiff__field}
            vtkdiff_config.getConfigParameter<std::string>("field");
        DBUG("vtkdiff will compare field '{:s}'.", field_name);

        std::vector<std::string> filenames;
        if (auto const regex_string =
                //! \ogs_file_param{prj__test_definition__vtkdiff__regex}
            vtkdiff_config.getConfigParameterOptional<std::string>("regex"))
        {
            // TODO: insert rank into regex for mpi case
            DBUG("vtkdiff regex is '{}'.", *regex_string);
            auto const regex = std::regex(*regex_string);
            for (auto const& p :
                 fs::directory_iterator(fs::path(reference_path)))
            {
                auto const filename = p.path().filename().string();
                if (std::regex_match(filename, regex))
                {
                    DBUG("        -> matched '{}'", filename);
                    filenames.push_back(filename);
                }
            }
        }
        else
        {
#ifdef USE_PETSC
            int rank;
            MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
            int mpi_size;
            MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
            std::string const& filename =
                MeshLib::IO::getVtuFileNameForPetscOutputWithoutExtension(
                    //! \ogs_file_param{prj__test_definition__vtkdiff__file}
                    vtkdiff_config.getConfigParameter<std::string>("file")) +
                "_" + std::to_string(rank) + ".vtu";
#else
            std::string const& filename =
                //! \ogs_file_param{prj__test_definition__vtkdiff__file}
                vtkdiff_config.getConfigParameter<std::string>("file");
#endif  // OGS_USE_PETSC
            filenames.push_back(filename);
        }

        if (empty(filenames))
        {
            OGS_FATAL(
                "No files from test definitions were added for tests but {} "
                "{:s} specified.",
                size(vtkdiff_configs),
                (size(vtkdiff_configs) == 1 ? "test was" : "tests were"));
        }

        auto const absolute_tolerance =
            //! \ogs_file_param{prj__test_definition__vtkdiff__absolute_tolerance}
            vtkdiff_config.getConfigParameter<std::string>("absolute_tolerance",
                                                           "");
        if (!absolute_tolerance.empty() &&
            !isConvertibleToDouble(absolute_tolerance))
        {
            OGS_FATAL(
                "The absolute tolerance value '{:s}' is not convertible to "
                "double.",
                absolute_tolerance);
        }
        std::string const absolute_tolerance_parameter =
            "--abs " + absolute_tolerance;
        auto const relative_tolerance =
            //! \ogs_file_param{prj__test_definition__vtkdiff__relative_tolerance}
            vtkdiff_config.getConfigParameter<std::string>("relative_tolerance",
                                                           "");
        if (!relative_tolerance.empty() &&
            !isConvertibleToDouble(relative_tolerance))
        {
            OGS_FATAL(
                "The relative tolerance value '{:s}' is not convertible to "
                "double.",
                relative_tolerance);
        }
        std::string const relative_tolerance_parameter =
            "--rel " + relative_tolerance;

        for (auto const& filename : filenames)
        {
            std::string const& output_filename =
                BaseLib::joinPaths(output_directory, filename);
            _output_files.push_back(output_filename);
            std::string const& reference_filename =
                BaseLib::joinPaths(reference_path, filename);

            //
            // Construct command line.
            //
            std::string command_line =
                vtkdiff + " -a " + safeString(field_name) + " -b " +
                safeString(field_name) + " " + safeString(reference_filename) +
                " " + safeString(output_filename) + " " +
                absolute_tolerance_parameter + " " +
                relative_tolerance_parameter;
            INFO("Will run '{:s}'", command_line);
            _command_lines.emplace_back(std::move(command_line));
        }
    }
}

bool TestDefinition::runTests() const
{
    std::vector<int> return_values;
    transform(begin(_command_lines), end(_command_lines),
              back_inserter(return_values),
              [](std::string const& command_line)
              {
                  int const return_value = std::system(command_line.c_str());
                  if (return_value != 0)
                  {
                      WARN("Return value {:d} was returned by '{:s}'.",
                           return_value, command_line);
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

std::size_t TestDefinition::numberOfTests() const
{
    return size(_command_lines);
}
}  // namespace ApplicationsLib
