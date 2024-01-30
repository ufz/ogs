/**
 * \file
 * \brief  Implementation of OpenGeoSys simulation application python module
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <spdlog/spdlog.h>
#include <tclap/CmdLine.h>

#include "../ogs.mesh/OGSMesh.h"
#include "Applications/ApplicationsLib/Simulation.h"
#include "Applications/ApplicationsLib/TestDefinition.h"
#include "BaseLib/DateTools.h"
#include "BaseLib/Error.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/RunTime.h"
#include "CommandLineArgumentParser.h"
#include "InfoLib/GitInfo.h"
#include "ogs_embedded_python.h"

std::unique_ptr<Simulation> simulation;

static constexpr int EXIT_ARGPARSE_FAILURE = 3;  // "mangled" TCLAP status
static constexpr int EXIT_ARGPARSE_EXIT_OK = 2;  // "mangled" TCLAP status
static_assert(EXIT_FAILURE == 1);
static_assert(EXIT_SUCCESS == 0);

int initOGS(std::vector<std::string>& argv_str)
{
    int argc = argv_str.size();
    char** argv = new char*[argc];
    for (int i = 0; i < argc; ++i)
    {
        argv[i] = argv_str[i].data();
    }

    CommandLineArguments cli_args;
    try
    {
        cli_args = parseCommandLineArguments(argc, argv, false);
    }
    catch (TCLAP::ArgException const& e)
    {
        ERR("Parsing the OGS commandline failed: {}", e.what());

        // "mangle" TCLAP's status
        return EXIT_ARGPARSE_FAILURE;
    }
    catch (TCLAP::ExitException const& e)
    {
        if (e.getExitStatus() == 0)
        {
            return EXIT_ARGPARSE_EXIT_OK;
        }

        // "mangle" TCLAP's status
        return EXIT_ARGPARSE_FAILURE;
    }

    BaseLib::initOGSLogger(cli_args.log_level);

    INFO("This is OpenGeoSys-6 version {:s}.",
         GitInfoLib::GitInfo::ogs_version);

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
        simulation = std::make_unique<Simulation>(argc, argv);
        simulation->initializeDataStructures(
            std::move(cli_args.project),
            std::move(cli_args.xml_patch_file_names),
            cli_args.reference_path_is_set, std::move(cli_args.reference_path),
            cli_args.nonfatal, std::move(cli_args.outdir),
            std::move(cli_args.mesh_dir), std::move(cli_args.script_dir),
            cli_args.write_prj);
    }
    catch (std::exception& e)
    {
        ERR("{}", e.what());
        ogs_status = EXIT_FAILURE;
    }

    INFO("OpenGeoSys is now initialized.");

    return ogs_status;
}

int executeSimulation()
{
    BaseLib::RunTime run_time;

    {
        auto const start_time = std::chrono::system_clock::now();
        auto const time_str = BaseLib::formatDate(start_time);
        INFO("OGS started on {:s}.", time_str);
    }

    auto ogs_status = EXIT_SUCCESS;

    try
    {
        run_time.start();
        bool solver_succeeded = simulation->executeSimulation();
        simulation->outputLastTimeStep();
        // TODO: test definition ?

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

    return ogs_status;
}

int executeTimeStep()
{
    auto ogs_status = EXIT_SUCCESS;
    try
    {
        bool solver_succeeded = simulation->executeTimeStep();
        ogs_status = solver_succeeded ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    catch (std::exception& e)
    {
        ERR("{}", e.what());
        ogs_status = EXIT_FAILURE;
    }
    return ogs_status;
}

double currentTime()
{
    return simulation->currentTime();
}

double endTime()
{
    return simulation->endTime();
}

OGSMesh getMesh(std::string const& name)
{
    return OGSMesh(simulation->getMesh(name));
}

void finalize()
{
    simulation.reset(nullptr);

    // TODO don't use global project directory, shared among different OGS
    // instances.
    // Unset project dir to make multiple OGS runs in one Python session
    // possible.
    BaseLib::unsetProjectDirectory();
}

/// To use this module import dependencies first:
///   import ogs.mesh as mesh
///   import ogs.simulator as sim
///
/// See also
/// https://github.com/pybind/pybind11/issues/1391#issuecomment-912642979
PYBIND11_MODULE(simulator, m)
{
    m.attr("__name__") = "ogs.simulator";
    m.doc() = "pybind11 ogs example plugin";
    m.def("initialize", &initOGS, "init OGS");
    m.def("currentTime", &currentTime, "get current OGS time");
    m.def("endTime", &endTime, "get end OGS time");
    m.def("executeSimulation", &executeSimulation, "execute OGS simulation");
    m.def("executeTimeStep", &executeTimeStep, "execute OGS time step");
    m.def("getMesh", &getMesh, "get unstructured grid from ogs");
    m.def("finalize", &finalize, "finalize OGS simulation");
}
