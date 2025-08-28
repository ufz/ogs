/**
 * \file
 * \brief  Implementation of OpenGeoSys simulation application python module
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <spdlog/spdlog.h>
#include <tclap/CmdLine.h>

#include <algorithm>

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

static constexpr int EXIT_ARGPARSE_FAILURE = 3;  // "mangled" TCLAP status
static constexpr int EXIT_ARGPARSE_EXIT_OK = 2;  // "mangled" TCLAP status
static_assert(EXIT_FAILURE == 1);
static_assert(EXIT_SUCCESS == 0);

// Needs to be exported, see
// https://pybind11.readthedocs.io/en/stable/advanced/misc.html#partitioning-code-over-multiple-extension-modules
class PYBIND11_EXPORT OGSSimulation
{
public:
    explicit OGSSimulation(std::vector<std::string>& argv_str)
    {
        INFO("OGSSimulation::OGSSimulation(std::vector<std::string>&)");

        int argc = argv_str.size();
        char** argv = new char*[argc];
        for (int i = 0; i < argc; ++i)
        {
            argv[i] = argv_str[i].data();
        }

        int cli_parse_status = EXIT_ARGPARSE_EXIT_OK;
        CommandLineArguments cli_args;
        try
        {
            cli_args = parseCommandLineArguments(argc, argv, false);
        }
        catch (TCLAP::ArgException const& e)
        {
            ERR("Parsing the OGS commandline failed: {}", e.what());

            // "mangle" TCLAP's status
            cli_parse_status = EXIT_ARGPARSE_FAILURE;
        }
        catch (TCLAP::ExitException const& e)
        {
            if (e.getExitStatus() == 0)
            {
                cli_parse_status = EXIT_ARGPARSE_EXIT_OK;
            }

            // "mangle" TCLAP's status
            cli_parse_status = EXIT_ARGPARSE_FAILURE;
        }

        BaseLib::initOGSLogger(cli_args.log_level);

        INFO(
            "This is OpenGeoSys-6 version {:s}. Log version: {:d}, Log level: "
            "{:s}.",
            GitInfoLib::GitInfo::ogs_version, 2, cli_args.log_level);

        BaseLib::createOutputDirectory(cli_args.outdir);

        {
            auto const start_time = std::chrono::system_clock::now();
            auto const time_str = BaseLib::formatDate(start_time);
            // todo ask Tobias: started vs starts
            INFO("OGS starts on {:s} in serial mode / Python embedded mode.",
                 time_str);
        }
        try
        {
            simulation = std::make_unique<Simulation>(argc, argv);
            simulation->initializeDataStructures(
                std::move(cli_args.project),
                std::move(cli_args.xml_patch_file_names),
                cli_args.reference_path_is_set,
                std::move(cli_args.reference_path), cli_args.nonfatal,
                std::move(cli_args.outdir), std::move(cli_args.mesh_dir),
                std::move(cli_args.script_dir), cli_args.write_prj);
        }
        catch (std::exception& e)
        {
            ERR("{}", e.what());
            ogs_status = EXIT_FAILURE;
        }
        INFO("OpenGeoSys is now initialized.");
    }

    int executeSimulation()
    {
        BaseLib::RunTime run_time;

        {
            auto const start_time = std::chrono::system_clock::now();
            auto const time_str = BaseLib::formatDate(start_time);
            INFO("OGS started on {:s} in serial mode.", time_str);
        }

        auto ogs_status = EXIT_SUCCESS;
        std::optional<ApplicationsLib::TestDefinition> test_definition{
            std::nullopt};

        try
        {
            run_time.start();
            bool solver_succeeded = simulation->executeSimulation();
            simulation->outputLastTimeStep();
            test_definition = simulation->getTestDefinition();

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

    double currentTime() { return simulation->currentTime(); }

    double endTime() { return simulation->endTime(); }

    OGSMesh& getMesh(std::string const& name)
    {
        auto const mesh_it = mesh_mapping.find(name);
        if (mesh_it != mesh_mapping.end())
        {
            INFO("found OGSMesh '{}' with address: {}", name,
                 fmt::ptr(&(mesh_it->second)));
            return mesh_it->second;
        }

        auto const& [it, success] =
            mesh_mapping.insert({name, OGSMesh(simulation->getMesh(name))});
        if (!success)
        {
            OGS_FATAL("Could not access mesh '{}'.", name);
        }
        INFO("insert OGSMesh '{}' with address: {}", name,
             fmt::ptr(&(it->second)));
        return it->second;
    }

    std::vector<std::string> getMeshNames()
    {
        return simulation->getMeshNames();
    }

    void finalize()
    {
        simulation->outputLastTimeStep();
        simulation.reset(nullptr);

        // Unset project dir to make multiple OGS runs in one Python session
        // possible.
        BaseLib::unsetProjectDirectory();
    }

private:
    int ogs_status = EXIT_SUCCESS;

    std::unique_ptr<Simulation> simulation;
    std::map<std::string, OGSMesh> mesh_mapping;
};

/// To use this module import dependencies first:
///   import ogs.mesh as mesh
///   import ogs.OGSSimulator as sim
///
/// See also
/// https://github.com/pybind/pybind11/issues/1391#issuecomment-912642979
PYBIND11_MODULE(OGSSimulator, m)
{
    m.attr("__name__") = "ogs.OGSSimulator";
    m.doc() = "pybind11 ogs plugin";

    pybind11::class_<OGSSimulation>(m, "OGSSimulation")
        .def(pybind11::init<std::vector<std::string>&>())
        .def("currentTime", &OGSSimulation::currentTime, "get current OGS time")
        .def("endTime", &OGSSimulation::endTime, "get end OGS time")
        .def("executeSimulation", &OGSSimulation::executeSimulation,
             "execute OGS simulation")
        .def("executeTimeStep", &OGSSimulation::executeTimeStep,
             "execute OGS time step")
        .def("getMesh", &OGSSimulation::getMesh,
             pybind11::return_value_policy::automatic_reference,
             pybind11::arg("name"), "get unstructured grid from ogs")
        .def("getMeshNames", &OGSSimulation::getMeshNames,
             "get names of all meshes from ogs")
        .def("finalize", &OGSSimulation::finalize, "finalize OGS simulation");
}
