/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include "CreateThermalTwoPhaseFlowWithPPProcess.h"

#include <cassert>

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
        curves)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "THERMAL_TWOPHASE_WITH_PP");

    DBUG("Create nonisothermal two-phase flow model.");
    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_THERMAL__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

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
    // diffusion coeff
    auto const& diff_coeff_b = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__TWOPHASE_FLOW_THERMAL__diffusion_coeff_component_b}
        "diffusion_coeff_component_b", parameters, 1, &mesh);
    auto const& diff_coeff_a = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__TWOPHASE_FLOW_THERMAL__diffusion_coeff_component_a}
        "diffusion_coeff_component_a", parameters, 1, &mesh);

    // Parameter for the density of the solid.

    auto const& density_solid = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__TWOPHASE_FLOW_THERMAL__density_solid}
        "density_solid", parameters, 1, &mesh);
    DBUG("Use '{:s}' as density_solid parameter.", density_solid.name);

    // Parameter for the latent heat of evaporation.
    auto const& latent_heat_evaporation = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__TWOPHASE_FLOW_THERMAL__latent_heat_evaporation}
        "latent_heat_evaporation", parameters, 1, &mesh);
    DBUG("Use '{:s}' as latent_heat_evaporation parameter.",
         latent_heat_evaporation.name);

    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_THERMAL__material_property}
    auto const& mat_config = config.getConfigSubtree("material_property");

    std::unique_ptr<ThermalTwoPhaseFlowWithPPMaterialProperties> material =
        createThermalTwoPhaseFlowWithPPMaterialProperties(
            mat_config, materialIDs(mesh), parameters);

    ThermalTwoPhaseFlowWithPPProcessData process_data{specific_body_force,
                                                      has_gravity,
                                                      mass_lumping,
                                                      diff_coeff_b,
                                                      diff_coeff_a,
                                                      density_solid,
                                                      latent_heat_evaporation,
                                                      std::move(material)};

    return std::make_unique<ThermalTwoPhaseFlowWithPPProcess>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables), mat_config,
        curves);
}

}  // namespace ThermalTwoPhaseFlowWithPP
}  // namespace ProcessLib
