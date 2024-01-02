/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include "CreateTwoPhaseFlowWithPrhoProcess.h"

#include <cassert>

#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/TwoPhaseFlowWithPrho/CreateTwoPhaseFlowPrhoMaterialProperties.h"
#include "ProcessLib/TwoPhaseFlowWithPrho/TwoPhaseFlowWithPrhoMaterialProperties.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "TwoPhaseFlowWithPrhoProcess.h"
#include "TwoPhaseFlowWithPrhoProcessData.h"

namespace ProcessLib
{
namespace TwoPhaseFlowWithPrho
{
void checkMPLProperties(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    std::array const required_medium_properties = {
        MaterialPropertyLib::permeability, MaterialPropertyLib::porosity};
    std::array const required_liquid_properties = {
        MaterialPropertyLib::viscosity, MaterialPropertyLib::density};
    std::array const required_gas_properties = {
        MaterialPropertyLib::viscosity, MaterialPropertyLib::density,
        MaterialPropertyLib::molar_mass};

    for (auto const& m : media)
    {
        checkRequiredProperties(*m.second, required_medium_properties);
        checkRequiredProperties(m.second->phase("AqueousLiquid"),
                                required_liquid_properties);
        checkRequiredProperties(m.second->phase("Gas"),
                                required_gas_properties);
    }
}

std::unique_ptr<Process> createTwoPhaseFlowWithPrhoProcess(
    std::string const& name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "TWOPHASE_FLOW_PRHO");

    DBUG("Create TwoPhaseFlowProcess with Prho model.");

    /// \section processvariablestpfwprho Process Variables

    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    /// Primary process variables as they appear in the global component vector:
    auto per_process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__TWOPHASE_FLOW_PRHO__process_variables__liquid_pressure}
         "liquid_pressure",
         //! \ogs_file_param_special{prj__processes__process__TWOPHASE_FLOW_PRHO__process_variables__overall_mass_density}
         "overall_mass_density"});
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    process_variables.push_back(std::move(per_process_variables));

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);
    /// \section parameterstpfwprho Process Parameters
    // Specific body force
    std::vector<double> const b =
        //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__specific_body_force}
        config.getConfigParameter<std::vector<double>>("specific_body_force");
    assert(!b.empty() && b.size() < 4);
    Eigen::VectorXd specific_body_force(b.size());
    bool const has_gravity = MathLib::toVector(b).norm() > 0;
    if (has_gravity)
    {
        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__mass_lumping}
    auto const mass_lumping = config.getConfigParameter<bool>("mass_lumping");
    // diffusion coeff
    auto const& diff_coeff_b = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__TWOPHASE_FLOW_PRHO__diffusion_coeff_component_b}
        "diffusion_coeff_component_b", parameters, 1, &mesh);
    auto const& diff_coeff_a = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__TWOPHASE_FLOW_PRHO__diffusion_coeff_component_a}
        "diffusion_coeff_component_a", parameters, 1, &mesh);
    auto const& temperature = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__TWOPHASE_FLOW_PRHO__temperature}
        "temperature", parameters, 1, &mesh);

    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_PRHO__material_property}
    auto const& mat_config = config.getConfigSubtree("material_property");

    auto const material_ids = materialIDs(mesh);
    if (material_ids != nullptr)
    {
        INFO("The twophase flow is in heterogeneous porous media.");
    }
    else
    {
        INFO("The twophase flow is in homogeneous porous media.");
    }

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);
    checkMPLProperties(media);

    std::unique_ptr<TwoPhaseFlowWithPrhoMaterialProperties> material =
        createTwoPhaseFlowPrhoMaterialProperties(mat_config, material_ids);

    TwoPhaseFlowWithPrhoProcessData process_data{
        specific_body_force, has_gravity,         mass_lumping,
        diff_coeff_b,        diff_coeff_a,        temperature,
        std::move(material), std::move(media_map)};

    return std::make_unique<TwoPhaseFlowWithPrhoProcess>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables));
}

}  // namespace TwoPhaseFlowWithPrho
}  // namespace ProcessLib
