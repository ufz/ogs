/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "MeshLib/Mesh.h"
#include "MeshLib/Vtk/VtkMappedMeshSource.h"

namespace InSituLib
{
vtkCPProcessor* Processor = NULL;

void Initialize(BaseLib::ConfigTree const & scripts_config, std::string const & path)
{
    if (Processor == NULL)
    {
        Processor = vtkCPProcessor::New();
        Processor->Initialize();
    }
    else
    {
        Processor->RemoveAllPipelines();
    }
    for (auto script_config : scripts_config.getConfigSubtreeList("script"))
    {
        auto scriptName = script_config.getConfigParameter<std::string>("name");
        INFO("Initializing in-situ script: %s", scriptName.c_str());
        std::stringstream ss;
        ss << path << scriptName;
        vtkNew<vtkCPPythonScriptPipeline> pipeline;
        pipeline->Initialize(ss.str().c_str());
        Processor->AddPipeline(pipeline.GetPointer());
    }
}
void Finalize()
{
    if (Processor)
    {
        Processor->Delete();
        Processor = NULL;
    }
}
void CoProcess(MeshLib::Mesh const& mesh, double const time,
               unsigned int const timeStep, bool const lastTimeStep)
{
    if (Processor == NULL)
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
        INFO("Start InSitu process.");
        vtkNew<MeshLib::VtkMappedMeshSource> vtkSource;
        vtkSource->SetMesh(&mesh);
        vtkSource->Update();
        dataDescription->GetInputDescriptionByName("input")->SetGrid(vtkSource->GetOutput());
        Processor->CoProcess(dataDescription.GetPointer());
        INFO("End InSitu process.");
    }
}
}  // namespace
