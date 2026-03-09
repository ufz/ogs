// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <spdlog/spdlog.h>
#include <tclap/CmdLine.h>

#include <algorithm>

#include "../ogs.OGSMesh/OGSMesh.h"
#include "Applications/ApplicationsLib/Simulation.h"
#include "Applications/ApplicationsLib/TestDefinition.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/DateTools.h"
#include "BaseLib/Error.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "BaseLib/RunTime.h"
#include "BaseLib/TCLAPArguments.h"
#include "CommandLineArgumentParser.h"
#include "InfoLib/GitInfo.h"

static constexpr int EXIT_ARGPARSE_FAILURE = 3;  // "mangled" TCLAP status
static constexpr int EXIT_ARGPARSE_EXIT_OK = 2;  // "mangled" TCLAP status
static_assert(EXIT_FAILURE == 1);
static_assert(EXIT_SUCCESS == 0);

#define OGS_ALWAYS_ASSERT(cond)                      \
    if (!(cond))                                     \
    {                                                \
        OGS_FATAL("OGS assertion failed {}", #cond); \
    }

std::pair<int, std::vector<char*>> toArgcArgv(
    std::vector<std::string>& argv_str)
{
    int argc = argv_str.size();
    std::vector<char*> argv_vec;
    argv_vec.reserve(argc + 1);
    for (auto& arg : argv_str)
    {
        argv_vec.push_back(arg.data());
    }
    argv_vec.push_back(nullptr);  // last entry must be a nullptr!

    return {argc, std::move(argv_vec)};
}

// Needs to be exported, see
// https://pybind11.readthedocs.io/en/stable/advanced/misc.html#partitioning-code-over-multiple-extension-modules
class PYBIND11_EXPORT OGSSimulation
{
public:
    explicit OGSSimulation(std::vector<std::string>& argv_str)
    {
        auto [argc, argv_vec] = toArgcArgv(argv_str);
        char** argv = argv_vec.data();

        mpi_setup.emplace(argc, argv);

        CommandLineArguments cli_args;
        try
        {
            cli_args = parseCommandLineArguments(argc, argv, false);
        }
        catch (TCLAP::ArgException const& e)
        {
            // TODO fragile interplay between (incomplete) simulation
            // initialization and OGS logger initialization
            BaseLib::initOGSLogger(BaseLib::defaultLogLevel());

            std::cerr << "Parsing the OGS commandline failed: " << e.what()
                      << '\n';

            // "mangle" TCLAP's status
            throw(e);
        }
        catch (TCLAP::ExitException const& e)
        {
            // TODO fragile interplay between (incomplete) simulation
            // initialization and OGS logger initialization
            BaseLib::initOGSLogger(BaseLib::defaultLogLevel());

            if (e.getExitStatus() == 0)
            {
                // --version/--help
                return;
            }

            throw(e);
        }

        BaseLib::initOGSLogger(cli_args.log_level);

        DBUG("OGSSimulation::OGSSimulation(std::vector<std::string>&)");

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
            simulation.reset();
            throw(e);
        }
        INFO("OpenGeoSys is now initialized.");
    }

    int executeSimulation()
    {
        OGS_ALWAYS_ASSERT(initialized());

        BaseLib::RunTime run_time;

        {
            auto const start_time = std::chrono::system_clock::now();
            auto const time_str = BaseLib::formatDate(start_time);
            INFO("OGS started on {:s} in serial mode.", time_str);
        }

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
        OGS_ALWAYS_ASSERT(initialized());

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

    double currentTime() const
    {
        OGS_ALWAYS_ASSERT(initialized());
        return simulation->currentTime();
    }

    double endTime() const
    {
        OGS_ALWAYS_ASSERT(initialized());
        return simulation->endTime();
    }

    OGSMesh& getMesh(std::string const& name)
    {
        OGS_ALWAYS_ASSERT(initialized());

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

    std::vector<std::string> getMeshNames() const
    {
        OGS_ALWAYS_ASSERT(initialized());

        return simulation->getMeshNames();
    }

    void finalize()
    {
        if (simulation)
        {
            simulation->outputLastTimeStep();
            simulation.reset(nullptr);
        }
        mpi_setup.reset();

        // Check for swallowed ConfigTree errors after Simulation destructor
        // runs. This catches configuration errors in objects destroyed at end
        // of scope.
        try
        {
            BaseLib::ConfigTree::assertNoSwallowedErrors();
        }
        catch (std::exception& e)
        {
            ERR("{}", e.what());
            throw;
        }
    }

    int status() const { return ogs_status; }

    bool initialized() const { return simulation != nullptr; }

private:
    int ogs_status = EXIT_SUCCESS;

    std::unique_ptr<Simulation> simulation;
    std::map<std::string, OGSMesh> mesh_mapping;
    std::optional<BaseLib::MPI::Setup> mpi_setup;
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
        .def("current_time", &OGSSimulation::currentTime,
             "get current OGS time")
        .def("end_time", &OGSSimulation::endTime, "get end OGS time")
        .def("execute_simulation", &OGSSimulation::executeSimulation,
             "execute OGS simulation")
        .def("execute_time_step", &OGSSimulation::executeTimeStep,
             "execute OGS time step")
        .def("mesh", &OGSSimulation::getMesh,
             pybind11::return_value_policy::automatic_reference,
             pybind11::arg("name"), "get unstructured grid from ogs")
        .def("mesh_names", &OGSSimulation::getMeshNames,
             "get names of all meshes from ogs")
        .def("close", &OGSSimulation::finalize, "finalize OGS simulation")
        .def_property_readonly("status", &OGSSimulation::status)
        .def_property_readonly(
            "initialized", &OGSSimulation::initialized,
            "Tells if the simulation object has been completely initialized.");
}
