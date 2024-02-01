/**
 * \file
 * \brief  Implementation of CommandLineArgumentParser.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CommandLineArgumentParser.h"

#include <tclap/CmdLine.h>

#include "BaseLib/FileTools.h"
#include "InfoLib/CMakeInfo.h"
#include "InfoLib/GitInfo.h"

CommandLineArguments parseCommandLineArguments(int argc, char* argv[],
                                               bool const exit_on_exception)
{
    // Parse CLI arguments.
    TCLAP::CmdLine cmd(
        "OpenGeoSys-6 software."
        "\n"
        "Copyright (c) 2012-2024, OpenGeoSys Community "
        "(https://www.opengeosys.org) "
        "\n"
        "Distributed under a Modified BSD License. "
        "See accompanying file LICENSE.txt or "
        "https://www.opengeosys.org/project/license"
        "\n"
        "version: " +
            GitInfoLib::GitInfo::ogs_version +
            "\n"
            "CMake arguments: " +
            CMakeInfoLib::CMakeInfo::cmake_args,
        ' ',
        GitInfoLib::GitInfo::ogs_version + "\n\n" +
            "CMake arguments: " + CMakeInfoLib::CMakeInfo::cmake_args);

    cmd.setExceptionHandling(exit_on_exception);

    TCLAP::ValueArg<std::string> log_level_arg(
        "l", "log-level",
        "the verbosity of logging messages: none, error, warn, info, "
        "debug, "
        "all",
        false,
#ifdef NDEBUG
        "info",
#else
        "all",
#endif
        "LOG_LEVEL");

#ifndef _WIN32  // TODO: On windows floating point exceptions are not handled
                // currently
    TCLAP::SwitchArg enable_fpe_arg("", "enable-fpe",
                                    "enables floating point exceptions");
#endif  // _WIN32
    TCLAP::SwitchArg unbuffered_cout_arg("", "unbuffered-std-out",
                                         "use unbuffered standard output");

    TCLAP::ValueArg<std::string> reference_path_arg(
        "r", "reference",
        "Run output result comparison after successful simulation "
        "comparing to all files in the given path. This requires test "
        "definitions to be present in the project file.",
        false, "", "PATH");

    TCLAP::UnlabeledValueArg<std::string> project_arg(
        "project-file",
        "Path to the ogs6 project file.",
        true,
        "",
        "PROJECT_FILE");

    TCLAP::MultiArg<std::string> xml_patch_files_arg(
        "p", "xml-patch",
        "the xml patch file(s) which is (are) applied (in the given order) "
        "to the PROJECT_FILE",
        false, "");

    TCLAP::ValueArg<std::string> outdir_arg("o", "output-directory",
                                            "the output directory to write to",
                                            false, "", "PATH");

    TCLAP::ValueArg<std::string> mesh_dir_arg(
        "m", "mesh-input-directory",
        "the directory where the meshes are read from", false, "", "PATH");

    TCLAP::ValueArg<std::string> script_dir_arg(
        "s", "script-input-directory",
        "the directory where script files (e.g. Python BCs) are read from",
        false, "", "PATH");

    TCLAP::SwitchArg write_prj_arg("",
                                   "write-prj",
                                   "Writes processed project file to output "
                                   "path / [prj_base_name]_processed.prj.");

    TCLAP::SwitchArg nonfatal_arg("",
                                  "config-warnings-nonfatal",
                                  "warnings from parsing the configuration "
                                  "file will not trigger program abortion");
    cmd.add(reference_path_arg);
    cmd.add(project_arg);
    cmd.add(xml_patch_files_arg);
    cmd.add(outdir_arg);
    cmd.add(mesh_dir_arg);
    cmd.add(script_dir_arg);
    cmd.add(write_prj_arg);
    cmd.add(log_level_arg);
    cmd.add(nonfatal_arg);
    cmd.add(unbuffered_cout_arg);
#ifndef _WIN32  // TODO: On windows floating point exceptions are not handled
                // currently
    cmd.add(enable_fpe_arg);
#endif  // _WIN32

    cmd.parse(argc, argv);

    CommandLineArguments cli_args;
    cli_args.reference_path = reference_path_arg.getValue();
    cli_args.reference_path_is_set = reference_path_arg.isSet();
    cli_args.project = project_arg.getValue();

    BaseLib::setProjectDirectory(BaseLib::extractPath(cli_args.project));

    cli_args.xml_patch_file_names = xml_patch_files_arg.getValue();
    cli_args.outdir = outdir_arg.getValue();
    cli_args.mesh_dir = mesh_dir_arg.getValue().empty()
                            ? BaseLib::getProjectDirectory()
                            : mesh_dir_arg.getValue();
    cli_args.script_dir = script_dir_arg.getValue().empty()
                              ? BaseLib::getProjectDirectory()
                              : script_dir_arg.getValue();
    cli_args.nonfatal = nonfatal_arg.getValue();
    cli_args.log_level = log_level_arg.getValue();
    cli_args.write_prj = write_prj_arg.getValue();

    // deactivate buffer for standard output if specified
    if (unbuffered_cout_arg.isSet())
    {
        std::cout.setf(std::ios::unitbuf);
    }
#ifndef _WIN32
    cli_args.enable_fpe_is_set = enable_fpe_arg.isSet();
#endif  // _WIN32

    return cli_args;
}
