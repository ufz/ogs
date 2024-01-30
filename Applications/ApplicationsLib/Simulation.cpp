/**
 * \file
 * \brief  Implementation of class Simulation
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Simulation.h"

#include <spdlog/spdlog.h>

#include "Applications/ApplicationsLib/LinearSolverLibrarySetup.h"
#include "Applications/ApplicationsLib/ProjectData.h"
#include "Applications/ApplicationsLib/TestDefinition.h"
#include "Applications/InSituLib/Adaptor.h"
#include "BaseLib/ConfigTreeUtil.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/PrjProcessing.h"
#include "MeshLib/Mesh.h"
#include "NumLib/NumericsConfig.h"
#include "ProcessLib/TimeLoop.h"

Simulation::Simulation(int argc, char* argv[])
    : linear_solver_library_setup(argc, argv),
#if defined(USE_PETSC)
      controller(vtkSmartPointer<vtkMPIController>::New()),
#endif
      test_definition{std::nullopt}
{
#if defined(USE_PETSC)
    controller->Initialize(&argc, &argv, 1);
    vtkMPIController::SetGlobalController(controller);

    {  // Can be called only after MPI_INIT.
        int mpi_rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &mpi_rank);
        spdlog::set_pattern(fmt::format("[{}] %^%l:%$ %v", mpi_rank));
    }
#endif
}

void Simulation::initializeDataStructures(
    std::string const& project,
    std::vector<std::string> const& xml_patch_file_names,
    bool const reference_path_is_set, std::string const& reference_path,
    bool const nonfatal, std::string const& outdir, std::string const& mesh_dir,
    std::string const& script_dir, bool const write_prj)
{
    INFO("Reading project file {}.",
         std::filesystem::absolute(project).string());

    std::stringstream prj_stream;
    BaseLib::prepareProjectFile(prj_stream, project, xml_patch_file_names,
                                write_prj, outdir);
    auto project_config = BaseLib::makeConfigTree(
        project, !nonfatal, "OpenGeoSysProject", prj_stream);

    if (!reference_path_is_set)
    {  // Ignore the test_definition section.
        project_config.ignoreConfigParameter("test_definition");
    }
    else
    {
        test_definition = ApplicationsLib::TestDefinition(
            //! \ogs_file_param{prj__test_definition}
            project_config.getConfigSubtree("test_definition"), reference_path,
            outdir);
        if (test_definition->numberOfTests() == 0)
        {
            OGS_FATAL(
                "No tests were constructed from the test definitions, "
                "but reference solutions path was given.");
        }

        INFO("Cleanup possible output files before running ogs.");
        BaseLib::removeFiles(test_definition->getOutputFiles());
    }
#ifdef OGS_USE_INSITU
    //! \ogs_file_param{prj__insitu}
    if (auto t = project_config.getConfigSubtreeOptional("insitu"))
    {
        InSituLib::Initialize(
            //! \ogs_file_param{prj__insitu__scripts}
            t->getConfigSubtree("scripts"),
            BaseLib::extractPath(project));
        isInsituConfigured = true;
    }
#else
    project_config.ignoreConfigParameter("insitu");
#endif

    project_data = std::make_unique<ProjectData>(project_config,
                                                 BaseLib::getProjectDirectory(),
                                                 outdir, mesh_dir, script_dir);

    INFO("Initialize processes.");
    for (auto& p : project_data->getProcesses())
    {
        p->initialize(project_data->getMedia());
    }

    // Check intermediately that config parsing went fine.
    checkAndInvalidate(project_config);
    BaseLib::ConfigTree::assertNoSwallowedErrors();

    auto& time_loop = project_data->getTimeLoop();
    time_loop.initialize();
}

double Simulation::currentTime() const
{
    auto const& time_loop = project_data->getTimeLoop();
    return time_loop.currentTime();
}

double Simulation::endTime() const
{
    auto const& time_loop = project_data->getTimeLoop();
    return time_loop.endTime();
}

bool Simulation::executeTimeStep()
{
    auto& time_loop = project_data->getTimeLoop();
    if (time_loop.currentTime() < time_loop.endTime())
    {
        auto const result = time_loop.executeTimeStep();
        if (time_loop.calculateNextTimeStep())
        {
            time_loop.outputLastTimeStep();
        }
        return result;
    }
    return false;
}

MeshLib::Mesh& Simulation::getMesh(std::string const& name)
{
    return project_data->getMesh(name);
}

bool Simulation::executeSimulation()
{
    INFO("Solve processes.");
    auto& time_loop = project_data->getTimeLoop();
    while (time_loop.currentTime() < time_loop.endTime())
    {
        time_loop.executeTimeStep();
        if (!time_loop.calculateNextTimeStep())
        {
            break;
        }
    }

    return time_loop.successful_time_step;
}

void Simulation::outputLastTimeStep() const
{
    auto const& time_loop = project_data->getTimeLoop();
    time_loop.outputLastTimeStep();
}

std::optional<ApplicationsLib::TestDefinition> Simulation::getTestDefinition()
    const
{
    return test_definition;
}

Simulation::~Simulation()
{
#ifdef OGS_USE_INSITU
    if (isInsituConfigured)
    {
        InSituLib::Finalize();
    }
#endif
#if defined(USE_PETSC)
    controller->Finalize(1);
#endif
}
