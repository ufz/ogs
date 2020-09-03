/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateHeatConductionProcess.h"

#include "HeatConductionProcess.h"
#include "HeatConductionProcessData.h"
#include "MaterialLib/MPL/CheckMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
namespace HeatConduction
{
void checkMPLProperties(
    MeshLib::Mesh const& mesh,
    MaterialPropertyLib::MaterialSpatialDistributionMap const& media_map)
{
    std::array const required_medium_properties = {
        MaterialPropertyLib::PropertyType::reference_temperature,
        MaterialPropertyLib::PropertyType::thermal_conductivity,
        MaterialPropertyLib::PropertyType::heat_capacity,
        MaterialPropertyLib::PropertyType::density};
    std::array<MaterialPropertyLib::PropertyType, 0> const empty{};

    MaterialPropertyLib::checkMaterialSpatialDistributionMap(
        mesh, media_map, required_medium_properties, empty, empty);
}

std::unique_ptr<Process> createHeatConductionProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "HEAT_CONDUCTION");

    DBUG("Create HeatConductionProcess.");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__HEAT_CONDUCTION__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    auto per_process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__HEAT_CONDUCTION__process_variables__process_variable}
         "process_variable"});
    process_variables.push_back(std::move(per_process_variables));

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    DBUG("Check the media properties of heat conduction process ...");
    checkMPLProperties(mesh, *media_map);
    DBUG("Media properties verified.");

    auto const mass_lumping =
        //! \ogs_file_param{prj__processes__process__HEAT_CONDUCTION__mass_lumping}
        config.getConfigParameter<bool>("mass_lumping", false);

    HeatConductionProcessData process_data{std::move(media_map), mass_lumping};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<HeatConductionProcess>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables));
}

}  // namespace HeatConduction
}  // namespace ProcessLib
