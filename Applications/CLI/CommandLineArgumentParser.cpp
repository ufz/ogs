/**
 * \brief  Implementation of CommandLineArgumentParser.
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <tclap/CmdLine.h>

#include "CommandLineArgumentParser.h"

#include "InfoLib/CMakeInfo.h"
#include "InfoLib/GitInfo.h"

CommandLineArgumentParser::CommandLineArgumentParser(int argc, char* argv[])
{
    // Parse CLI arguments.
    TCLAP::CmdLine cmd(
        "OpenGeoSys-6 software.\n"
        "Copyright (c) 2012-2022, OpenGeoSys Community "
        "(http://www.opengeosys.org) "
        "Distributed under a Modified BSD License. "
        "See accompanying file LICENSE.txt or "
        "http://www.opengeosys.org/project/license\n"
        "version: " +
            GitInfoLib::GitInfo::ogs_version + "\n" +
            "CMake arguments: " + CMakeInfoLib::CMakeInfo::cmake_args,
        ' ',
        GitInfoLib::GitInfo::ogs_version + "\n\n" +
            "CMake arguments: " + CMakeInfoLib::CMakeInfo::cmake_args);

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

    TCLAP::SwitchArg nonfatal_arg("",
                                  "config-warnings-nonfatal",
                                  "warnings from parsing the configuration "
                                  "file will not trigger program abortion");
    cmd.add(reference_path_arg);
    cmd.add(project_arg);
    cmd.add(xml_patch_files_arg);
    cmd.add(outdir_arg);
    cmd.add(log_level_arg);
    cmd.add(nonfatal_arg);
    cmd.add(unbuffered_cout_arg);
#ifndef _WIN32  // TODO: On windows floating point exceptions are not handled
                // currently
    cmd.add(enable_fpe_arg);
#endif  // _WIN32

    cmd.parse(argc, argv);

    reference_path = reference_path_arg.getValue();
    reference_path_is_set = reference_path_arg.isSet();
    project = project_arg.getValue();
    xml_patch_file_names = xml_patch_files_arg.getValue();
    outdir = outdir_arg.getValue();
    nonfatal = nonfatal_arg.getValue();
    log_level = log_level_arg.getValue();

    // deactivate buffer for standard output if specified
    if (unbuffered_cout_arg.isSet())
    {
        std::cout.setf(std::ios::unitbuf);
    }
#ifndef _WIN32
    enable_fpe_is_set = enable_fpe_arg.isSet();
#endif  // _WIN32
}
