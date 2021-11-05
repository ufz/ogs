/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateSteadyStateDiffusion.h"

#include "BaseLib/FileTools.h"
#include "MaterialLib/MPL/CheckMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "SteadyStateDiffusion.h"
#include "SteadyStateDiffusionData.h"

namespace ProcessLib
{
namespace SteadyStateDiffusion
{
void checkMPLProperties(
    MeshLib::Mesh const& mesh,
    MaterialPropertyLib::MaterialSpatialDistributionMap const& media_map)
{
    std::array const required_medium_properties = {
        MaterialPropertyLib::PropertyType::reference_temperature,
        MaterialPropertyLib::PropertyType::diffusion};
    std::array<MaterialPropertyLib::PropertyType, 0> const empty{};

    MaterialPropertyLib::checkMaterialSpatialDistributionMap(
        mesh, media_map, required_medium_properties, empty, empty);
}

std::unique_ptr<Process> createSteadyStateDiffusion(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "STEADY_STATE_DIFFUSION");

    DBUG("Create SteadyStateDiffusion.");

    /// \section processvariablesssd Process Variables

    //! \ogs_file_param{prj__processes__process__STEADY_STATE_DIFFUSION__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    /// Primary process variables as they appear in the global component vector:
    auto per_process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__STEADY_STATE_DIFFUSION__process_variables__process_variable}
         "process_variable"});
    process_variables.push_back(std::move(per_process_variables));

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    DBUG("Check the media properties of steady state diffusion process ...");
    checkMPLProperties(mesh, *media_map);
    DBUG("Media properties verified.");

    SteadyStateDiffusionData process_data{std::move(media_map)};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    /// \section parametersssd Process Parameters
    std::unique_ptr<ProcessLib::SurfaceFluxData> surfaceflux;
    auto calculatesurfaceflux_config =
        //! \ogs_file_param{prj__processes__process__calculatesurfaceflux}
        config.getConfigSubtreeOptional("calculatesurfaceflux");
    if (calculatesurfaceflux_config)
    {
        surfaceflux = ProcessLib::SurfaceFluxData::createSurfaceFluxData(
            *calculatesurfaceflux_config, meshes);
    }

    return std::make_unique<SteadyStateDiffusion>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        std::move(surfaceflux));
}

}  // namespace SteadyStateDiffusion
}  // namespace ProcessLib
