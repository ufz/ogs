/**
 * \file
 * \date   2014-08-04
 * \brief  Implementation of OpenGeoSys simulation application
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <spdlog/spdlog.h>
#include <tclap/CmdLine.h>

#include <chrono>
#include <sstream>

#include "CommandLineArgumentParser.h"

#ifndef _WIN32
#ifdef __APPLE__
#ifdef __SSE__
#include <xmmintrin.h>
#endif  // __SSE__
#else
#include <cfenv>
#endif  // __APPLE__
#endif  // _WIN32

#include "Applications/ApplicationsLib/Simulation.h"
#include "Applications/ApplicationsLib/TestDefinition.h"
#include "BaseLib/DateTools.h"
#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"
#include "BaseLib/RunTime.h"
#include "InfoLib/GitInfo.h"

#ifndef _WIN32  // On windows this command line option is not present.
void enableFloatingPointExceptions()
{
#ifdef __APPLE__
#ifdef __SSE__
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
#endif  // __SSE__
#else
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif  // __APPLE__
}
#endif  // _WIN32

int main(int argc, char* argv[])
{
    CommandLineArguments cli_arg = parseCommandLineArguments(argc, argv);
    BaseLib::initOGSLogger(cli_arg.log_level);
#ifndef _WIN32  // TODO: On windows floating point exceptions are not handled
    if (cli_arg.enable_fpe_is_set)
    {
        enableFloatingPointExceptions();
    }
#endif  // _WIN32

    INFO("This is OpenGeoSys-6 version {:s}.",
         GitInfoLib::GitInfo::ogs_version);

    BaseLib::RunTime run_time;

    {
        auto const start_time = std::chrono::system_clock::now();
        auto const time_str = BaseLib::formatDate(start_time);
        INFO("OGS started on {:s}.", time_str);
    }

    std::optional<ApplicationsLib::TestDefinition> test_definition{
        std::nullopt};
    auto ogs_status = EXIT_SUCCESS;

    try
    {
        Simulation simulation(argc, argv);
        run_time.start();
        simulation.initializeDataStructures(
            std::move(cli_arg.project), std::move(cli_arg.xml_patch_file_names),
            cli_arg.reference_path_is_set, std::move(cli_arg.reference_path),
            cli_arg.nonfatal, std::move(cli_arg.outdir),
            std::move(cli_arg.mesh_dir), std::move(cli_arg.script_dir),
            cli_arg.write_prj);
        bool solver_succeeded = simulation.executeSimulation();
        simulation.outputLastTimeStep();
        test_definition = simulation.getTestDefinition();

        INFO("[time] Execution took {:g} s.", run_time.elapsed());
        ogs_status = solver_succeeded ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    catch (std::exception& e)
    {
        ERR("{}", e.what());
        ogs_status = EXIT_FAILURE;
    }

    {
        auto const end_time = std::chrono::system_clock::now();
        auto const time_str = BaseLib::formatDate(end_time);
        INFO("OGS terminated on {:s}.", time_str);
    }

    if (ogs_status == EXIT_FAILURE)
    {
        ERR("OGS terminated with error.");
        return EXIT_FAILURE;
    }

    if (!test_definition)
    {
        // There are no tests, so just exit;
        return ogs_status;
    }

    INFO("");
    INFO("##########################################");
    INFO("# Running tests                          #");
    INFO("##########################################");
    INFO("");
    if (!test_definition->runTests())
    {
        ERR("One of the tests failed.");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
