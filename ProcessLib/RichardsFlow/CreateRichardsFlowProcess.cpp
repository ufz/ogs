/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateRichardsFlowProcess.h"

#include "MaterialLib/MPL/CheckMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "ParameterLib/ConstantParameter.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "RichardsFlowProcess.h"
#include "RichardsFlowProcessData.h"

namespace ProcessLib
{
namespace RichardsFlow
{
void checkMPLProperties(
    MeshLib::Mesh const& mesh,
    MaterialPropertyLib::MaterialSpatialDistributionMap const& media_map)
{
    std::array const required_medium_properties = {
        MaterialPropertyLib::PropertyType::reference_temperature,
        MaterialPropertyLib::PropertyType::saturation,
        MaterialPropertyLib::PropertyType::relative_permeability};

    std::array const required_liquid_properties = {
        MaterialPropertyLib::PropertyType::viscosity,
        MaterialPropertyLib::PropertyType::density};

    std::array const required_solid_properties = {
        MaterialPropertyLib::PropertyType::porosity,
        MaterialPropertyLib::PropertyType::storage,
        MaterialPropertyLib::PropertyType::permeability};

    MaterialPropertyLib::checkMaterialSpatialDistributionMap(
        mesh, media_map, required_medium_properties, required_solid_properties,
        required_liquid_properties);
}

std::unique_ptr<Process> createRichardsFlowProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        /*curves*/,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "RICHARDS_FLOW");

    DBUG("Create RichardsFlowProcess.");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__RICHARDS_FLOW__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    auto per_process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__RICHARDS_FLOW__process_variables__process_variable}
         "process_variable"});
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    process_variables.push_back(std::move(per_process_variables));

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    // Specific body force
    std::vector<double> const b =
        //! \ogs_file_param{prj__processes__process__RICHARDS_FLOW__specific_body_force}
        config.getConfigParameter<std::vector<double>>("specific_body_force");
    assert(!b.empty() && b.size() < 4);
    Eigen::VectorXd specific_body_force(b.size());
    bool const has_gravity = MathLib::toVector(b).norm() > 0;
    if (has_gravity)
    {
        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    // has mass lumping
    //! \ogs_file_param{prj__processes__process__RICHARDS_FLOW__mass_lumping}
    auto mass_lumping = config.getConfigParameter<bool>("mass_lumping");

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    DBUG("Check the media properties of RichardsFlow process ...");
    checkMPLProperties(mesh, *media_map);
    DBUG("Media properties verified.");

    RichardsFlowProcessData process_data{
        std::move(media_map), specific_body_force, has_gravity, mass_lumping};

    return std::make_unique<RichardsFlowProcess>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables));
}

}  // namespace RichardsFlow
}  // namespace ProcessLib
