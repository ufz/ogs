/**
 * \date   2014-08-04
 * \brief  Implementation of OpenGeoSys simulation application
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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

#ifdef USE_PETSC
#include <vtkMPIController.h>
#include <vtkSmartPointer.h>
#endif

#include "Applications/ApplicationsLib/Simulation.h"
#include "Applications/ApplicationsLib/TestDefinition.h"
#include "BaseLib/DateTools.h"
#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"
#include "BaseLib/RunTime.h"
#include "InfoLib/GitInfo.h"

#ifdef OGS_USE_PYTHON
#include "ogs_embedded_python.h"
#endif

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

void setConsoleLogLevel(std::string const& log_level)
{
    BaseLib::setConsoleLogLevel(log_level);
    spdlog::set_pattern("%^%l:%$ %v");
    spdlog::set_error_handler(
        [](const std::string& msg)
        {
            std::cerr << "spdlog error: " << msg << std::endl;
            OGS_FATAL("spdlog logger error occurred.");
        });
}

int main(int argc, char* argv[])
{
    CommandLineArgumentParser cli_arg(argc, argv);
    setConsoleLogLevel(cli_arg.log_level);
#ifndef _WIN32  // TODO: On windows floating point exceptions are not handled
    if (cli_arg.enable_fpe_is_set)
    {
        enableFloatingPointExceptions();
    }
#endif  // _WIN32

    INFO("This is OpenGeoSys-6 version {:s}.",
         GitInfoLib::GitInfo::ogs_version);

#ifdef OGS_USE_PYTHON
    pybind11::scoped_interpreter guard = ApplicationsLib::setupEmbeddedPython();
    (void)guard;
#endif

    BaseLib::RunTime run_time;

    {
        auto const start_time = std::chrono::system_clock::now();
        auto const time_str = BaseLib::formatDate(start_time);
        INFO("OGS started on {:s}.", time_str);
    }

    std::unique_ptr<ApplicationsLib::TestDefinition> test_definition;
    auto ogs_status = EXIT_SUCCESS;

    try
    {
        Simulation ogs(argc, argv);
        run_time.start();
        ogs.initializeDataStructures(
            std::move(cli_arg.project), std::move(cli_arg.xml_patch_file_names),
            cli_arg.reference_path_is_set, std::move(cli_arg.reference_path),
            cli_arg.nonfatal, std::move(cli_arg.outdir));
        bool solver_succeeded = ogs.executeSimulation();

        INFO("[time] Execution took {:g} s.", run_time.elapsed());
        ogs_status = solver_succeeded ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    catch (std::exception& e)
    {
        ERR(e.what());
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

    if (test_definition == nullptr)
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
