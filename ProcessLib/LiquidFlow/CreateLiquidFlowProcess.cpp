/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 * Created on August 19, 2016, 1:30 PM
 */
#include "CreateLiquidFlowProcess.h"

#include <algorithm>

#include "MaterialLib/MPL/CheckMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/PhysicalConstant.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "LiquidFlowProcess.h"

namespace ProcessLib
{
namespace LiquidFlow
{
void checkMPLProperties(
    MeshLib::Mesh const& mesh,
    MaterialPropertyLib::MaterialSpatialDistributionMap const& media_map)
{
    std::array const required_medium_properties = {
        MaterialPropertyLib::reference_temperature,
        MaterialPropertyLib::PropertyType::permeability,
        MaterialPropertyLib::PropertyType::porosity,
        MaterialPropertyLib::PropertyType::storage};

    std::array const required_liquid_properties = {
        MaterialPropertyLib::PropertyType::viscosity,
        MaterialPropertyLib::PropertyType::density};

    std::array<MaterialPropertyLib::PropertyType, 0> const
        required_solid_properties{};

    MaterialPropertyLib::checkMaterialSpatialDistributionMap(
        mesh, media_map, required_medium_properties, required_solid_properties,
        required_liquid_properties);
}

std::unique_ptr<Process> createLiquidFlowProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
    std::string const& output_directory,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "LIQUID_FLOW");

    DBUG("Create LiquidFlowProcess.");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    auto per_process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__LIQUID_FLOW__process_variables__process_variable}
         "process_variable"});
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    process_variables.push_back(std::move(per_process_variables));

    SecondaryVariableCollection secondary_variables;


    ProcessLib::createSecondaryVariables(config, secondary_variables);

    // Get the gravity vector for the Darcy velocity
    //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__darcy_gravity}
    auto const& darcy_g_config = config.getConfigSubtree("darcy_gravity");
    const auto gravity_axis_id_input =
        //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__darcy_gravity__axis_id}
        darcy_g_config.getConfigParameter<int>("axis_id");
    if (gravity_axis_id_input >= static_cast<int>(mesh.getDimension()))
    {
        OGS_FATAL(
            "The gravity axis must be a number between 0 and one less than the "
            "mesh dimension, which is {:d}. Read gravity axis {:d} from input "
            "file.",
            mesh.getDimension(), gravity_axis_id_input);
    }
    const auto g =
        //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__darcy_gravity__g}
        darcy_g_config.getConfigParameter<double>("g");
    if (g < 0.)
    {
        OGS_FATAL("Gravity magnitude must be non-negative.");
    }
    const int gravity_axis_id = (g == 0.) ? -1 : gravity_axis_id_input;

    std::unique_ptr<ProcessLib::SurfaceFluxData> surfaceflux;
    auto calculatesurfaceflux_config =
        //! \ogs_file_param{prj__processes__process__calculatesurfaceflux}
        config.getConfigSubtreeOptional("calculatesurfaceflux");
    if (calculatesurfaceflux_config)
    {
        surfaceflux = ProcessLib::SurfaceFluxData::createSurfaceFluxData(
            *calculatesurfaceflux_config, meshes);
    }

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    DBUG("Check the media properties of LiquidFlow process ...");
    checkMPLProperties(mesh, *media_map);
    DBUG("Media properties verified.");

    LiquidFlowData process_data{std::move(media_map), gravity_axis_id, g};

    return std::make_unique<LiquidFlowProcess>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        std::move(surfaceflux));
}

}  // namespace LiquidFlow
}  // namespace ProcessLib
