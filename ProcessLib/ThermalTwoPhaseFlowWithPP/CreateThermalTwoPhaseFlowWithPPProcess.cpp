/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include "CreateThermalTwoPhaseFlowWithPPProcess.h"

#include <cassert>

#include "BaseLib/Functional.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "ProcessLib/Parameter/ConstantParameter.h"
#include "ProcessLib/ThermalTwoPhaseFlowWithPP/CreateThermalTwoPhaseFlowWithPPMaterialProperties.h"
#include "ProcessLib/ThermalTwoPhaseFlowWithPP/ThermalTwoPhaseFlowWithPPMaterialProperties.h"
#include "ProcessLib/Utils/ParseSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "ThermalTwoPhaseFlowWithPPProcess.h"
#include "ThermalTwoPhaseFlowWithPPProcessData.h"

namespace ProcessLib
{
namespace ThermalTwoPhaseFlowWithPP
{
std::unique_ptr<Process> createThermalTwoPhaseFlowWithPPProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
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

    auto process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__TWOPHASE_FLOW_THERMAL__process_variables__gas_pressure}
         "gas_pressure",
         //! \ogs_file_param_special{prj__processes__process__TWOPHASE_FLOW_THERMAL__process_variables__capillary_pressure}
         "capillary_pressure",
         //! \ogs_file_param_special{prj__processes__process__TWOPHASE_FLOW_THERMAL__process_variables__temperature}
         "temperature"});

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"TwoPhaseFlow_pressure"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);
    // Specific body force
    std::vector<double> const b =
        //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_THERMAL__specific_body_force}
        config.getConfigParameter<std::vector<double>>("specific_body_force");
    assert(b.size() > 0 && b.size() < 4);
    Eigen::VectorXd specific_body_force(b.size());
    bool const has_gravity = MathLib::toVector(b).norm() > 0;
    if (has_gravity)
        std::copy_n(b.data(), b.size(), specific_body_force.data());

    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_THERMAL__mass_lumping}
    auto mass_lumping = config.getConfigParameter<bool>("mass_lumping");
    // diffusion coeff
    auto& diff_coeff_b = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__TWOPHASE_FLOW_THERMAL__diffusion_coeff_component_b}
        "diffusion_coeff_component_b", parameters, 1);
    auto& diff_coeff_a = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__TWOPHASE_FLOW_THERMAL__diffusion_coeff_component_a}
        "diffusion_coeff_component_a", parameters, 1);

    // Parameter for the density of the solid.

    auto& density_solid = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__TWOPHASE_FLOW_THERMAL__density_solid}
        "density_solid", parameters, 1);
    DBUG("Use \'%s\' as density_solid parameter.", density_solid.name.c_str());

    // Parameter for the latent heat of evaporation.
    auto& latent_heat_evaporation = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__TWOPHASE_FLOW_THERMAL__latent_heat_evaporation}
        "latent_heat_evaporation", parameters, 1);
    DBUG("Use \'%s\' as latent_heat_evaporation parameter.",
         latent_heat_evaporation.name.c_str());

    //! \ogs_file_param{prj__processes__process__TWOPHASE_FLOW_THERMAL__material_property}
    auto const& mat_config = config.getConfigSubtree("material_property");

    std::unique_ptr<ThermalTwoPhaseFlowWithPPMaterialProperties>
        material = nullptr;

    if (mesh.getProperties().existsPropertyVector<int>("MaterialIDs"))
    {
        INFO("The twophase flow is in heterogeneous porous media.");
        auto const& mat_ids =
            mesh.getProperties().getPropertyVector<int>("MaterialIDs");
        material = createThermalTwoPhaseFlowWithPPMaterialProperties(mat_config,
                                                                     *mat_ids);
    }
    else
    {
        INFO("The twophase flow is in homogeneous porous media.");
        MeshLib::Properties dummy_property;
        auto const& dummy_property_vector =
            dummy_property.createNewPropertyVector<int>(
                "MaterialIDs", MeshLib::MeshItemType::Cell, 1);
        material = createThermalTwoPhaseFlowWithPPMaterialProperties(
            mat_config, *dummy_property_vector);
    }

    ThermalTwoPhaseFlowWithPPProcessData process_data{specific_body_force,
                                                      has_gravity,
                                                      mass_lumping,
                                                      diff_coeff_b,
                                                      diff_coeff_a,
                                                      density_solid,
                                                      latent_heat_evaporation,
                                                      std::move(material)};

    return std::make_unique<ThermalTwoPhaseFlowWithPPProcess>(
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller),
        mat_config, curves);
}

}  // end of namespace
}  // end of namespace
