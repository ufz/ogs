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
#include <ctime>
#include <sstream>

// ThirdParty/tclap
#include "tclap/CmdLine.h"

// BaseLib
#include "BaseLib/BuildInfo.h"
#include "BaseLib/ConfigTreeUtil.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/RunTime.h"

#include "Applications/ApplicationsLib/LinearSolverLibrarySetup.h"
#include "Applications/ApplicationsLib/LogogSetup.h"
#include "Applications/ApplicationsLib/ProjectData.h"
#include "Applications/ApplicationsLib/UncoupledProcessesTimeLoop.h"

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

    cmd.parse(argc, argv);

    ApplicationsLib::LogogSetup logog_setup;
    logog_setup.setLevel(log_level_arg.getValue());

    BaseLib::RunTime run_time;
    run_time.start();

    {
        auto const start_time_sys = std::chrono::system_clock::to_time_t(
            std::chrono::system_clock::now());
        std::ostringstream sstr;
        sstr << std::put_time(std::localtime(&start_time_sys), "%F %T %z");
        INFO("OGS started on %s.", sstr.str().c_str());
    }

    std::chrono::steady_clock::now();

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

            // Create processes.
            project.buildProcesses();

            BaseLib::ConfigTree::assertNoSwallowedErrors();

            INFO("Initialize processes.");
            for (auto p_it = project.processesBegin();
                 p_it != project.processesEnd(); ++p_it)
            {
                (*p_it)->initialize();
            }

            BaseLib::ConfigTree::assertNoSwallowedErrors();

            INFO("Solve processes.");

            auto& time_loop = project.getTimeLoop();
            solver_succeeded = time_loop.loop(project);
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
        auto const end_time_sys = std::chrono::system_clock::to_time_t(
            std::chrono::system_clock::now());
        std::ostringstream sstr;
        sstr << std::put_time(std::localtime(&end_time_sys), "%F %T %z");
        INFO("OGS terminated on %s.", sstr.str().c_str());
        INFO("[time] Execution took %g s.", run_time.elapsed());
    }

    return ogs_status;
}
