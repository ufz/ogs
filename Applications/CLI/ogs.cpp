/**
 * \date   2014-08-04
 * \brief  Implementation of OpenGeoSys simulation application
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <chrono>
#include <tclap/CmdLine.h>

// BaseLib
#include "BaseLib/BuildInfo.h"
#include "BaseLib/ConfigTreeUtil.h"
#include "BaseLib/DateTools.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/RunTime.h"

#include "Applications/ApplicationsLib/LinearSolverLibrarySetup.h"
#include "Applications/ApplicationsLib/LogogSetup.h"
#include "Applications/ApplicationsLib/ProjectData.h"
#include "ProcessLib/UncoupledProcessesTimeLoop.h"

#include "NumLib/NumericsConfig.h"


int main(int argc, char *argv[])
{
    // Parse CLI arguments.
    TCLAP::CmdLine cmd("OpenGeoSys-6 software.\n"
            "Copyright (c) 2012-2016, OpenGeoSys Community "
            "(http://www.opengeosys.org) "
            "Distributed under a Modified BSD License. "
            "See accompanying file LICENSE.txt or "
            "http://www.opengeosys.org/project/license\n"
            "version: " + BaseLib::BuildInfo::git_describe,
        ' ',
        BaseLib::BuildInfo::git_describe);

    TCLAP::UnlabeledValueArg<std::string> project_arg(
        "project-file",
        "Path to the ogs6 project file.",
        true,
        "",
        "PROJECT FILE");
    cmd.add(project_arg);

    TCLAP::ValueArg<std::string> outdir_arg(
        "o", "output-directory",
        "the output directory to write to",
        false,
        "",
        "output directory");
    cmd.add(outdir_arg);

    TCLAP::ValueArg<std::string> log_level_arg(
        "l", "log-level",
        "the verbosity of logging messages: none, error, warn, info, debug, all",
        false,
        "all",
        "log level");
    cmd.add(log_level_arg);

    TCLAP::SwitchArg nonfatal_arg("",
        "config-warnings-nonfatal",
        "warnings from parsing the configuration file will not trigger program abortion");
    cmd.add(nonfatal_arg);

    TCLAP::SwitchArg unbuffered_cout_arg("",
        "unbuffered-std-out",
        "use unbuffered standard output");
    cmd.add(unbuffered_cout_arg);

    cmd.parse(argc, argv);

    // deactivate buffer for standard output if specified
    if (unbuffered_cout_arg.isSet())
        std::cout.setf(std::ios::unitbuf);

    ApplicationsLib::LogogSetup logog_setup;
    logog_setup.setLevel(log_level_arg.getValue());

    INFO("This is OpenGeoSys-6 version %s.",
         BaseLib::BuildInfo::git_describe.c_str());

    BaseLib::RunTime run_time;
    run_time.start();

    {
        auto const start_time = std::chrono::system_clock::now();
        auto const time_str = BaseLib::formatDate(start_time);
        INFO("OGS started on %s.", time_str.c_str());
    }

    auto ogs_status = EXIT_SUCCESS;

    try
    {
        bool solver_succeeded = false;
        {
            ApplicationsLib::LinearSolverLibrarySetup
                linear_solver_library_setup(argc, argv);

            auto project_config = BaseLib::makeConfigTree(
                project_arg.getValue(), !nonfatal_arg.getValue(),
                "OpenGeoSysProject");

            ProjectData project(*project_config,
                                BaseLib::extractPath(project_arg.getValue()),
                                outdir_arg.getValue());

            // Check intermediately that config parsing went fine.
            project_config.checkAndInvalidate();
            BaseLib::ConfigTree::assertNoSwallowedErrors();

            BaseLib::ConfigTree::assertNoSwallowedErrors();

            INFO("Initialize processes.");
            for (auto& p : project.getProcesses())
            {
                p.second->initialize();
            }

            BaseLib::ConfigTree::assertNoSwallowedErrors();

            INFO("Solve processes.");

            auto& time_loop = project.getTimeLoop();
            solver_succeeded = time_loop.loop();
        }  // This nested scope ensures that everything that could possibly
           // possess a ConfigTree is destructed before the final check below is
           // done.

        BaseLib::ConfigTree::assertNoSwallowedErrors();

        ogs_status = solver_succeeded ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (std::exception& e) {
        ERR(e.what());
        ogs_status = EXIT_FAILURE;
    }

    {
        auto const end_time = std::chrono::system_clock::now();
        auto const time_str = BaseLib::formatDate(end_time);
        INFO("OGS terminated on %s.", time_str.c_str());
    }

    INFO("[time] Execution took %g s.", run_time.elapsed());

    return ogs_status;
}
