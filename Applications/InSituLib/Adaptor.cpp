/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Adaptor.h"

#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>

#include <filesystem>

#include "BaseLib/ConfigTree.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Vtk/VtkMappedMeshSource.h"

namespace InSituLib
{
vtkCPProcessor* Processor = nullptr;

void Initialize(BaseLib::ConfigTree const& scripts_config,
                std::string const& path)
{
    if (Processor == nullptr)
    {
        Processor = vtkCPProcessor::New();
        Processor->Initialize();
    }
    else
    {
        Processor->RemoveAllPipelines();
    }
    //! \ogs_file_param{prj__insitu__scripts__script}
    for (auto script_config : scripts_config.getConfigSubtreeList("script"))
    {
        auto scriptPath = std::filesystem::path(
            //! \ogs_file_param{prj__insitu__scripts__script__name}
            script_config.getConfigParameter<std::string>("name"));
        if (scriptPath.is_relative())
        {
            scriptPath = std::filesystem::path(path) / scriptPath;
        }
        if (!std::filesystem::exists(scriptPath))
        {
            ERR("In-situ script {:s} does not exist!", scriptPath.string());
        }
        INFO("Initializing in-situ script: {:s}", scriptPath.string());
        vtkNew<vtkCPPythonScriptPipeline> pipeline;
        pipeline->Initialize(scriptPath.c_str());
        Processor->AddPipeline(pipeline.GetPointer());
    }
}
void Finalize()
{
    if (Processor)
    {
        Processor->Delete();
        Processor = nullptr;
    }
}
void CoProcess(MeshLib::Mesh const& mesh, double const time,
               unsigned int const timeStep, bool const lastTimeStep,
               std::string output_directory)
{
    if (Processor == nullptr)
        return;

    vtkNew<vtkCPDataDescription> dataDescription;
    dataDescription->AddInput("input");
    dataDescription->SetTimeData(time, timeStep);
    if (lastTimeStep == true)
    {
        // assume that we want to all the pipelines to execute if it
        // is the last time step.
        dataDescription->ForceOutputOn();
    }
    if (Processor->RequestDataDescription(dataDescription.GetPointer()) != 0)
    {
        INFO("Start InSitu process: timestep #{:d} (t={:g}, last={:d})",
             timeStep, time, lastTimeStep);
        vtkNew<MeshLib::VtkMappedMeshSource> vtkSource;
        vtkSource->SetMesh(&mesh);
        vtkSource->Update();
        dataDescription->GetInputDescriptionByName("input")->SetGrid(
            vtkSource->GetOutput());
        auto const cwd = std::filesystem::current_path();
        std::filesystem::current_path(std::move(output_directory));
        Processor->CoProcess(dataDescription.GetPointer());
        std::filesystem::current_path(cwd);
        INFO("End InSitu process.");
    }
}
}  // namespace InSituLib
