// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "TestDefinition.h"

#include <spdlog/fmt/fmt.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/MPI.h"

#ifdef USE_PETSC
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
        ERR("The given string '{:s}' is not convertible to double.", s);
        return false;
    }
    if (pos != s.size())
    {
        ERR("Only {:d} characters were used for double conversion of string "
            "'{:s}'",
            pos, s);
        return false;
    }

    if (std::isnan(value))
    {
        ERR("The given string '{:s}' results in a NaN value.", s);
        return false;
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

/// Tries to find a diff tool executable by testing 'path/tool --version' calls
/// for various paths.
std::string findDiffTool(std::string const& executable_name,
                         std::string const& environment_variable_name)
{
    // Try to read the environment variable.
    if (const char* diff_tool_exe_environment_variable =
            std::getenv(environment_variable_name.c_str()))
    {
        std::string const diff_tool_exe{diff_tool_exe_environment_variable};
        DBUG("{:s} set to {:s}.", environment_variable_name, diff_tool_exe);

        //
        // Sanity checks.
        //
        {  // Check the base name.
            auto const& base_name =
                BaseLib::extractBaseNameWithoutExtension(diff_tool_exe);
            if (base_name != executable_name)
            {
                OGS_FATAL(
                    "The {:s} environment variable does not point to '{:s}'. "
                    "{:s}='{:s}'",
                    environment_variable_name, executable_name,
                    environment_variable_name, diff_tool_exe);
            }
        }
        {  // Diff tool must exist.
            if (!BaseLib::IsFileExisting(diff_tool_exe))
            {
                OGS_FATAL("The {:s} points to a non-existing file. {:s}='{:s}'",
                          environment_variable_name, environment_variable_name,
                          diff_tool_exe);
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
            std::system((diff_tool_exe + " --version").c_str());
        if (return_value == 0)
        {
            return diff_tool_exe;
        }
        WARN(
            "Calling {:s} from the {:s} environment variable didn't work as "
            "expected. Return value was {:d}.",
            diff_tool_exe, environment_variable_name, return_value);
    }

    std::vector<std::string> const paths = {"", "bin"};
    auto const path =
        find_if(begin(paths), end(paths),
                [&executable_name](std::string const& path)
                {
                    int const return_value =
                        // TODO (naumov) replace system call with output
                        // consuming call as in an above todo comment.
                        std::system((BaseLib::joinPaths(path, executable_name) +
                                     " --version")
                                        .c_str());
                    return return_value == 0;
                });
    if (path == end(paths))
    {
        OGS_FATAL("{:s} not found.", executable_name);
    }
    return BaseLib::joinPaths(*path, executable_name);
}

std::string formatPathForDiffTool(std::string const& filename)
{
#if _WIN32
    // VTK does not handle Windows long paths:
    // https://gitlab.kitware.com/vtk/vtk/-/blob/master/Utilities/KWSys/vtksys/SystemTools.cxx#L1519-1521
    // Workaround is to make the path absolute and prefix with a special
    // marker and put everything in quotes, see
    // https://learn.microsoft.com/en-us/windows/win32/fileio/maximum-file-path-limitation?tabs=powershell
    auto const& long_path_indicator = R"(\\?\)";

    auto const absolute_filename = std::filesystem::absolute(filename).string();
    return fmt::format("\"{}{}\"", long_path_indicator, absolute_filename);
#else
    return safeString(filename);
#endif
}

void checkTolerance(std::string const& tolerance_name,
                    std::string const& tolerance)
{
    if (!tolerance.empty() && !isConvertibleToDouble(tolerance))
    {
        OGS_FATAL(
            "The {:s} tolerance value '{:s}' is not convertible to double.",
            tolerance_name, tolerance);
    }
}

/// Pairs of reference and output file names to be compared.
using FilePairs = std::vector<std::pair<std::string, std::string>>;

/// Collects all file names in reference_path matching the regex.
FilePairs findFilesMatchingRegex(std::string const& regex_string,
                                 std::string const& reference_path)
{
    DBUG("regex is '{}'.", regex_string);
    FilePairs filenames;
    auto const regex = std::regex(regex_string);
    for (auto const& p : std::filesystem::directory_iterator(
             std::filesystem::path(reference_path)))
    {
        auto const filename = p.path().filename().string();
        if (std::regex_match(filename, regex))
        {
            DBUG("        -> matched '{}'", filename);
            filenames.emplace_back(filename, filename);
        }
    }
    return filenames;
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

    // Construct command lines for each entry.
    //! \ogs_file_param{prj__test_definition__vtkdiff}
    auto const& vtkdiff_configs = config_tree.getConfigSubtreeList("vtkdiff");
    //! \ogs_file_param{prj__test_definition__xdmfdiff}
    auto const& xdmfdiff_configs = config_tree.getConfigSubtreeList("xdmfdiff");
    _command_lines.reserve(vtkdiff_configs.size() + xdmfdiff_configs.size());

    // Constructs one command line per file name pair from the entries common
    // to all diff tools and appends it together with the corresponding output
    // file.
    auto const append_tests =
        [&](std::string const& tool, std::string const& tool_path,
            std::size_t const number_of_configs, FilePairs const& filenames,
            std::string const& field_name,
            std::string const& reference_field_name,
            std::string const& absolute_tolerance,
            std::string const& relative_tolerance,
            std::string const& extra_options)
    {
        if (filenames.empty())
        {
            OGS_FATAL(
                "No files from test definitions were added for tests but {} "
                "{:s} specified.",
                number_of_configs,
                (number_of_configs == 1 ? "test was" : "tests were"));
        }
        checkTolerance("absolute", absolute_tolerance);
        checkTolerance("relative", relative_tolerance);

        for (auto const& [reference_file, output_file] : filenames)
        {
            auto const output_filename =
                BaseLib::joinPaths(output_directory, output_file);
            _output_files.push_back(output_filename);

            std::string command_line = fmt::format(
                "{} -a {} -b {} {} {} --abs {} --rel {}{}", tool_path,
                safeString(reference_field_name), safeString(field_name),
                formatPathForDiffTool(
                    BaseLib::joinPaths(reference_path, reference_file)),
                formatPathForDiffTool(output_filename), absolute_tolerance,
                relative_tolerance, extra_options);
            INFO("Will run '{:s}'", command_line);
            _command_lines.push_back({tool, std::move(command_line)});
        }
    };

    std::string const vtkdiff =
        vtkdiff_configs.empty() ? "" : findDiffTool("vtkdiff", "VTKDIFF_EXE");
    for (auto const& vtkdiff_config : vtkdiff_configs)
    {
        std::string const& field_name =
            //! \ogs_file_param{prj__test_definition__vtkdiff__field}
            vtkdiff_config.getConfigParameter<std::string>("field");
        DBUG("vtkdiff will compare field '{:s}'.", field_name);
        std::string const reference_field_name =
            vtkdiff_config
                //! \ogs_file_param{prj__test_definition__vtkdiff__reference_field}
                .getConfigParameterOptional<std::string>("reference_field")
                .value_or(field_name);

        FilePairs filenames;
        if (auto const regex_string =
                //! \ogs_file_param{prj__test_definition__vtkdiff__regex}
            vtkdiff_config.getConfigParameterOptional<std::string>("regex"))
        {
            // TODO: insert rank into regex for mpi case
            filenames = findFilesMatchingRegex(*regex_string, reference_path);
        }
        else
        {
            std::string filename =
                //! \ogs_file_param{prj__test_definition__vtkdiff__file}
                vtkdiff_config.getConfigParameter<std::string>("file");
            std::string reference_filename =
                vtkdiff_config
                    //! \ogs_file_param{prj__test_definition__vtkdiff__reference_file}
                    .getConfigParameterOptional<std::string>("reference_file")
                    .value_or(filename);
#ifdef USE_PETSC
            BaseLib::MPI::Mpi mpi;
            if (mpi.size > 1)
            {
                filename =
                    MeshLib::IO::getVtuFileNameForPetscOutputWithoutExtension(
                        filename) +
                    "_" + std::to_string(mpi.rank) + ".vtu";
                reference_filename =
                    MeshLib::IO::getVtuFileNameForPetscOutputWithoutExtension(
                        reference_filename) +
                    "_" + std::to_string(mpi.rank) + ".vtu";
            }
#endif  // OGS_USE_PETSC
            filenames.emplace_back(reference_filename, filename);
        }

        auto const absolute_tolerance =
            //! \ogs_file_param{prj__test_definition__vtkdiff__absolute_tolerance}
            vtkdiff_config.getConfigParameter<std::string>("absolute_tolerance",
                                                           "");
        auto const relative_tolerance =
            //! \ogs_file_param{prj__test_definition__vtkdiff__relative_tolerance}
            vtkdiff_config.getConfigParameter<std::string>("relative_tolerance",
                                                           "");

        append_tests("vtkdiff", vtkdiff, vtkdiff_configs.size(), filenames,
                     field_name, reference_field_name, absolute_tolerance,
                     relative_tolerance, "");
    }

    std::string const xdmfdiff = xdmfdiff_configs.empty()
                                     ? ""
                                     : findDiffTool("xdmfdiff", "XDMFDIFF_EXE");
    for (auto const& xdmfdiff_config : xdmfdiff_configs)
    {
        std::string const& field_name =
            //! \ogs_file_param{prj__test_definition__xdmfdiff__field}
            xdmfdiff_config.getConfigParameter<std::string>("field");
        DBUG("xdmfdiff will compare field '{:s}'.", field_name);
        std::string const reference_field_name =
            xdmfdiff_config
                //! \ogs_file_param{prj__test_definition__xdmfdiff__reference_field}
                .getConfigParameterOptional<std::string>("reference_field")
                .value_or(field_name);

        FilePairs filenames;
        if (auto const regex_string =
                //! \ogs_file_param{prj__test_definition__xdmfdiff__regex}
            xdmfdiff_config.getConfigParameterOptional<std::string>("regex"))
        {
            filenames = findFilesMatchingRegex(*regex_string, reference_path);
        }
        else
        {
            std::string const filename =
                //! \ogs_file_param{prj__test_definition__xdmfdiff__file}
                xdmfdiff_config.getConfigParameter<std::string>("file");
            std::string const reference_filename =
                xdmfdiff_config
                    //! \ogs_file_param{prj__test_definition__xdmfdiff__reference_file}
                    .getConfigParameterOptional<std::string>("reference_file")
                    .value_or(filename);
            filenames.emplace_back(reference_filename, filename);
        }

        auto const absolute_tolerance =
            //! \ogs_file_param{prj__test_definition__xdmfdiff__absolute_tolerance}
            xdmfdiff_config.getConfigParameter<std::string>(
                "absolute_tolerance", "");
        auto const relative_tolerance =
            //! \ogs_file_param{prj__test_definition__xdmfdiff__relative_tolerance}
            xdmfdiff_config.getConfigParameter<std::string>(
                "relative_tolerance", "");

        std::string const timestep =
            //! \ogs_file_param{prj__test_definition__xdmfdiff__timestep}
            xdmfdiff_config.getConfigParameter<std::string>("timestep");
        std::string const reference_timestep =
            xdmfdiff_config
                //! \ogs_file_param{prj__test_definition__xdmfdiff__reference_timestep}
                .getConfigParameterOptional<std::string>("reference_timestep")
                .value_or(timestep);

        append_tests("xdmfdiff", xdmfdiff, xdmfdiff_configs.size(), filenames,
                     field_name, reference_field_name, absolute_tolerance,
                     relative_tolerance,
                     fmt::format(" --timestep-a {} --timestep-b {}",
                                 reference_timestep, timestep));
    }
}

bool TestDefinition::runTests() const
{
    return runCommandLines(_command_lines, true);
}

bool TestDefinition::runTests(std::string_view const diff_tool_name) const
{
    std::vector<CommandLine> command_lines;
    copy_if(begin(_command_lines), end(_command_lines),
            back_inserter(command_lines),
            [diff_tool_name](CommandLine const& command_line)
            { return command_line.diff_tool_name == diff_tool_name; });

    return runCommandLines(command_lines, true);
}

bool TestDefinition::runTestsExcluding(
    std::string_view const diff_tool_name) const
{
    std::vector<CommandLine> command_lines;
    copy_if(begin(_command_lines), end(_command_lines),
            back_inserter(command_lines),
            [diff_tool_name](CommandLine const& command_line)
            { return command_line.diff_tool_name != diff_tool_name; });

    return runCommandLines(command_lines, false);
}

bool TestDefinition::hasTests(std::string_view const diff_tool_name) const
{
    return any_of(begin(_command_lines), end(_command_lines),
                  [diff_tool_name](CommandLine const& command_line)
                  { return command_line.diff_tool_name == diff_tool_name; });
}

bool TestDefinition::runCommandLines(
    std::vector<CommandLine> const& command_lines, bool const fail_on_empty)
{
    std::vector<int> return_values;
    transform(begin(command_lines), end(command_lines),
              back_inserter(return_values),
              [](CommandLine const& command_line)
              {
                  INFO("---------- {:s} begin ----------",
                       command_line.diff_tool_name);
                  int const return_value =
                      std::system(command_line.command_line.c_str());
                  if (return_value != 0)
                  {
                      WARN("Value {:d} was returned by '{:s}'.", return_value,
                           command_line.command_line);
                  }
                  INFO("---------- {:s} end ----------\n",
                       command_line.diff_tool_name);
                  return return_value;
              });
    return (!return_values.empty() || !fail_on_empty) &&
           all_of(begin(return_values), end(return_values),
                  [](int const return_value) { return return_value == 0; });
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
