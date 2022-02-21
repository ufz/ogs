/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include "CreateThermalTwoPhaseFlowWithPPProcess.h"

#include <cassert>

#include "MaterialLib/MPL/CheckMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "ParameterLib/ConstantParameter.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/ThermalTwoPhaseFlowWithPP/CreateThermalTwoPhaseFlowWithPPMaterialProperties.h"
#include "ProcessLib/ThermalTwoPhaseFlowWithPP/ThermalTwoPhaseFlowWithPPMaterialProperties.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "ThermalTwoPhaseFlowWithPPProcess.h"
#include "ThermalTwoPhaseFlowWithPPProcessData.h"

namespace ProcessLib
{
namespace ThermalTwoPhaseFlowWithPP
{
void checkMPLProperties(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    std::array const required_property_medium = {
        MaterialPropertyLib::PropertyType::porosity,
        MaterialPropertyLib::PropertyType::permeability,
        MaterialPropertyLib::PropertyType::saturation};

    std::array const required_property_solid_phase = {
        MaterialPropertyLib::PropertyType::specific_heat_capacity,
        MaterialPropertyLib::PropertyType::density};

    std::array const required_property_liquid_phase = {
        MaterialPropertyLib::PropertyType::viscosity,
        MaterialPropertyLib::PropertyType::specific_heat_capacity,
        MaterialPropertyLib::PropertyType::density};

    std::array const required_property_gas_phase = {
        MaterialPropertyLib::PropertyType::viscosity};

    std::array const required_property_vapor_component = {
        MaterialPropertyLib::specific_heat_capacity};

    std::array const required_property_vapor_component = {
        MaterialPropertyLib::diffusion};

    std::array const required_property_vapor_component = {
        MaterialPropertyLib::specific_latent_heat};

    std::array const required_property_dry_air_component = {
        MaterialPropertyLib::specific_heat_capacity};

    for (auto const& m : media)
    {
        auto const& gas_phase = m.second->phase("Gas");
        checkRequiredProperties(*m.second, required_property_medium);
        checkRequiredProperties(gas_phase, required_property_gas_phase);
        checkRequiredProperties(m.second->phase("AqueousLiquid"),
                                required_property_liquid_phase);
        checkRequiredProperties(m.second->phase("Solid"),
                                required_property_solid_phase);

        // TODO (BM): should use index to identify components (same for impl.h)
        checkRequiredProperties(gas_phase.component("w"),
                                required_property_vapor_component);
        checkRequiredProperties(gas_phase.component("a"),
                                required_property_dry_air_component);
    }
}

std::unique_ptr<Process> createThermalTwoPhaseFlowWithPPProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "THERMAL_TWOPHASE_WITH_PP");

    /// \section processvariables Process Variables
    DBUG("Create nonisothermal two-phase flow model.");
    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_THERMAL__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    /// Primary process variables as they appear in the global component vector:
    auto per_process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__TWOPHASE_FLOW_THERMAL__process_variables__gas_pressure}
         "gas_pressure",
         //! \ogs_file_param_special{prj__processes__process__TWOPHASE_FLOW_THERMAL__process_variables__capillary_pressure}
         "capillary_pressure",
         //! \ogs_file_param_special{prj__processes__process__TWOPHASE_FLOW_THERMAL__process_variables__temperature}
         "temperature"});
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    process_variables.push_back(std::move(per_process_variables));

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);
    /// \section parameters Process Parameters
    // Specific body force
    std::vector<double> const b =
        //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_THERMAL__specific_body_force}
        config.getConfigParameter<std::vector<double>>("specific_body_force");
    assert(!b.empty() && b.size() < 4);
    Eigen::VectorXd specific_body_force(b.size());
    bool const has_gravity = MathLib::toVector(b).norm() > 0;
    if (has_gravity)
    {
        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_THERMAL__mass_lumping}
    auto mass_lumping = config.getConfigParameter<bool>("mass_lumping");

    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_THERMAL__material_property}
    auto const& mat_config = config.getConfigSubtree("material_property");

    std::unique_ptr<ThermalTwoPhaseFlowWithPPMaterialProperties> material =
        createThermalTwoPhaseFlowWithPPMaterialProperties(
            mat_config, materialIDs(mesh), parameters);

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    DBUG(
        "Check the media properties of ThermalTwoPhaseFlowWithPP  process ...");
    checkMPLProperties(media);
    DBUG("Media properties verified.");

    ThermalTwoPhaseFlowWithPPProcessData process_data{std::move(media_map),
                                                      specific_body_force,
                                                      has_gravity,
                                                      mass_lumping,
                                                      std::move(material)};

    return std::make_unique<ThermalTwoPhaseFlowWithPPProcess>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables), mat_config,
        curves);
}

}  // namespace ThermalTwoPhaseFlowWithPP
}  // namespace ProcessLib
