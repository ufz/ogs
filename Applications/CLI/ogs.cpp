// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <pybind11/pybind11.h>
#include <spdlog/spdlog.h>
#include <tclap/CmdLine.h>

#include <algorithm>
#include <chrono>
#include <csignal>
#include <iostream>
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
#include "BaseLib/FileTools.h"
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

#include <spdlog/sinks/null_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "BaseLib/MPI.h"

void signalHandler(int signum)
{
    auto const end_time = std::chrono::system_clock::now();
    auto const time_str = BaseLib::formatDate(end_time);

    ERR("Simulation aborted on {:s}. Received signal: {:d}.", time_str, signum);
    exit(signum);
}

void initializeLogger([[maybe_unused]] bool const all_ranks_log)
{
#if defined(USE_PETSC)
    int mpi_rank;
    MPI_Comm_rank(BaseLib::MPI::OGS_COMM_WORLD, &mpi_rank);
    int world_size;
    MPI_Comm_size(BaseLib::MPI::OGS_COMM_WORLD, &world_size);

    if (all_ranks_log)
    {
        if (world_size > 1)
        {
            spdlog::set_pattern(fmt::format("[{}] %^%l:%$ %v", mpi_rank));
        }
        // else untouched
    }
    else  // only rank 0 logs
    {
        // set_pattern is untouched
        spdlog::drop_all();
        if (mpi_rank > 0)
        {
            BaseLib::console = spdlog::create<spdlog::sinks::null_sink_st>(
                "ogs");  // do not log
        }
        else  // rank 0
        {
            auto console = BaseLib::console =
                spdlog::stdout_color_st("ogs");  // st for performance
        }
    }

    {
        auto const start_time = std::chrono::system_clock::now();
        auto const time_str = BaseLib::formatDate(start_time);
        INFO("OGS started on {:s} with MPI. MPI processes: {:d}.", time_str,
             world_size);
    }

#else  // defined(USE_PETSC)

    {
        auto const start_time = std::chrono::system_clock::now();
        auto const time_str = BaseLib::formatDate(start_time);
        INFO("OGS started on {:s} in serial mode.", time_str);
    }

#endif
}

int main(int argc, char* argv[])
{
    CommandLineArguments cli_arg = parseCommandLineArguments(argc, argv);

    // Initialize MPI
    // also in python hook
    // check tools
    BaseLib::MPI::Setup mpi_setup(argc, argv);
    BaseLib::initOGSLogger(cli_arg.log_level);
    initializeLogger(cli_arg.log_parallel);

    signal(SIGINT, signalHandler);   // CTRL+C
    signal(SIGTERM, signalHandler);  // pkill -SIGTERM <process_id> , It is NOT
                                     // possible to catch SIGKILL

#ifndef _WIN32  // TODO: On windows floating point exceptions are not handled
    if (cli_arg.enable_fpe_is_set)
    {
        enableFloatingPointExceptions();
    }
#endif  // _WIN32

    INFO(
        "This is OpenGeoSys-6 version {:s}. Log version: {:d}, Log level: "
        "{:s}.",
        GitInfoLib::GitInfo::ogs_version, 2, cli_arg.log_level);
    BaseLib::createOutputDirectory(cli_arg.outdir);

    std::optional<ApplicationsLib::TestDefinition> test_definition{
        std::nullopt};
    auto ogs_status = EXIT_SUCCESS;

    try
    {
        Simulation simulation(argc, argv);

        BaseLib::RunTime run_time;
        run_time.start();

        bool solver_succeeded = false;
        try
        {
            simulation.initializeDataStructures(
                std::move(cli_arg.project),
                std::move(cli_arg.xml_patch_file_names),
                cli_arg.reference_path_is_set,
                std::move(cli_arg.reference_path), cli_arg.nonfatal,
                std::move(cli_arg.outdir), std::move(cli_arg.mesh_dir),
                std::move(cli_arg.script_dir), cli_arg.write_prj);
            solver_succeeded = simulation.executeSimulation();
            simulation.outputLastTimeStep();
            test_definition = simulation.getTestDefinition();
        }
        catch (pybind11::error_already_set const& e)
        {
            OGS_FATAL("Python exception thrown: {}", e.what());
        }
        if (solver_succeeded)
        {
            INFO("[time] Simulation completed. It took {:g} s.",
                 run_time.elapsed());
        }
        else
        {
            INFO("[time] Simulation failed. It took {:g} s.",
                 run_time.elapsed());
        }

        ogs_status = solver_succeeded ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    catch (std::exception& e)
    {
        ERR("{}", e.what());
        ogs_status = EXIT_FAILURE;
    }

    if (ogs_status == EXIT_FAILURE)
    {
        auto const end_time = std::chrono::system_clock::now();
        auto const time_str = BaseLib::formatDate(end_time);
        ERR("OGS terminated with error on {:s}.", time_str);
        return EXIT_FAILURE;
    }

    return Simulation::runTestDefinitions(test_definition);
}
