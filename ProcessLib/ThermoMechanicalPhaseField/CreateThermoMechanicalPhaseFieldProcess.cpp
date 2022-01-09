/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateThermoMechanicalPhaseFieldProcess.h"

#include <cassert>

#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "ThermoMechanicalPhaseFieldProcess.h"
#include "ThermoMechanicalPhaseFieldProcessData.h"

namespace ProcessLib
{
namespace ThermoMechanicalPhaseField
{
template <int DisplacementDim>
std::unique_ptr<Process> createThermoMechanicalPhaseFieldProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "THERMO_MECHANICAL_PHASE_FIELD");
    DBUG("Create ThermoMechanicalPhaseFieldProcess.");

    INFO(
        "Solve the coupling with the staggered scheme,"
        "which is the only option for TM-Phasefield in the current code");

    /// \section processvariablestmpf Process Variables

    //! \ogs_file_param{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    int heat_conduction_process_id = 0;
    int mechanics_related_process_id = 1;
    int phase_field_process_id = 2;

    /// Primary process variables as they appear in the global component vector:
    auto process_variable_T = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__process_variables__temperature}
         "temperature"});
    process_variables.push_back(std::move(process_variable_T));
    ProcessVariable* variable_T =
        &process_variables[process_variables.size() - 1][0].get();

    auto process_variable_u = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__process_variables__displacement}
         "displacement"});
    process_variables.push_back(std::move(process_variable_u));
    ProcessVariable* variable_u =
        &process_variables[process_variables.size() - 1][0].get();
    auto process_variable_ph = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__process_variables__phasefield}
         "phasefield"});
    process_variables.push_back(std::move(process_variable_ph));
    ProcessVariable* variable_ph =
        &process_variables[process_variables.size() - 1][0].get();

    DBUG("Associate displacement with process variable '{:s}'.",
         variable_u->getName());

    if (variable_u->getNumberOfGlobalComponents() != DisplacementDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '{:s}' is different "
            "from the displacement dimension: got {:d}, expected {:d}",
            variable_u->getName(),
            variable_u->getNumberOfGlobalComponents(),
            DisplacementDim);
    }

    DBUG("Associate phase field with process variable '{:s}'.",
         variable_ph->getName());
    if (variable_ph->getNumberOfGlobalComponents() != 1)
    {
        OGS_FATAL(
            "Phasefield process variable '{:s}' is not a scalar variable but "
            "has {:d} components.",
            variable_ph->getName(),
            variable_ph->getNumberOfGlobalComponents());
    }

    DBUG("Associate temperature with process variable '{:s}'.",
         variable_T->getName());
    if (variable_T->getNumberOfGlobalComponents() != 1)
    {
        OGS_FATAL(
            "Temperature process variable '{:s}' is not a scalar variable but "
            "has {:d} components.",
            variable_T->getName(),
            variable_T->getNumberOfGlobalComponents());
    }

    /// \section parameterstmpf Process Parameters
    auto solid_constitutive_relations =
        MaterialLib::Solids::createConstitutiveRelations<DisplacementDim>(
            parameters, local_coordinate_system, config);

    auto const phasefield_parameters_config =
        //! \ogs_file_param{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__phasefield_parameters}
        config.getConfigSubtree("phasefield_parameters");

    auto const thermal_parameters_config =
        //! \ogs_file_param{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__thermal_parameters}
        config.getConfigSubtree("thermal_parameters");

    // Residual stiffness
    auto const& residual_stiffness = ParameterLib::findParameter<double>(
        phasefield_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__phasefield_parameters__residual_stiffness}
        "residual_stiffness", parameters, 1, &mesh);
    DBUG("Use '{:s}' as residual stiffness.", residual_stiffness.name);

    // Crack resistance
    auto const& crack_resistance = ParameterLib::findParameter<double>(
        phasefield_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__phasefield_parameters__crack_resistance}
        "crack_resistance", parameters, 1, &mesh);
    DBUG("Use '{:s}' as crack resistance.", crack_resistance.name);

    // Crack length scale
    auto const& crack_length_scale = ParameterLib::findParameter<double>(
        phasefield_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__phasefield_parameters__crack_length_scale}
        "crack_length_scale", parameters, 1, &mesh);
    DBUG("Use '{:s}' as crack length scale.", crack_length_scale.name);

    // Kinetic coefficient
    auto const& kinetic_coefficient = ParameterLib::findParameter<double>(
        phasefield_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__phasefield_parameters__kinetic_coefficient}
        "kinetic_coefficient", parameters, 1, &mesh);
    DBUG("Use '{:s}' as kinetic coefficient.", kinetic_coefficient.name);

    // Solid density
    auto const& solid_density = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__reference_solid_density}
        "solid_density", parameters, 1, &mesh);
    DBUG("Use '{:s}' as solid density parameter.", solid_density.name);

    // Linear thermal expansion coefficient
    auto const& linear_thermal_expansion_coefficient =
        ParameterLib::findParameter<double>(
            thermal_parameters_config,
            //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__thermal_parameters__linear_thermal_expansion_coefficient}
            "linear_thermal_expansion_coefficient", parameters, 1, &mesh);
    DBUG("Use '{:s}' as linear thermal expansion coefficient.",
         linear_thermal_expansion_coefficient.name);

    // Specific heat capacity
    auto const& specific_heat_capacity = ParameterLib::findParameter<double>(
        thermal_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__thermal_parameters__specific_heat_capacity}
        "specific_heat_capacity", parameters, 1, &mesh);
    DBUG("Use '{:s}' as specific heat capacity.", specific_heat_capacity.name);

    // Thermal conductivity
    auto const& thermal_conductivity = ParameterLib::findParameter<double>(
        thermal_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__thermal_parameters__thermal_conductivity}
        "thermal_conductivity", parameters, 1, &mesh);
    DBUG("Use '{:s}' as thermal conductivity parameter.",
         thermal_conductivity.name);
    // Residual thermal conductivity
    auto const& residual_thermal_conductivity = ParameterLib::findParameter<
        double>(
        thermal_parameters_config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__thermal_parameters__residual_thermal_conductivity}
        "residual_thermal_conductivity", parameters, 1, &mesh);
    DBUG("Use '{:s}' as residual thermal conductivity parameter.",
         residual_thermal_conductivity.name);
    // Reference temperature
    const auto reference_temperature =
        //! \ogs_file_param{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__reference_temperature}
        config.getConfigParameter<double>("reference_temperature");

    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__specific_body_force}
            config.getConfigParameter<std::vector<double>>(
                "specific_body_force");
        if (specific_body_force.size() != DisplacementDim)
        {
            OGS_FATAL(
                "The size of the specific body force vector does not match the "
                "displacement dimension. Vector size is {:d}, displacement "
                "dimension is {:d}",
                specific_body_force.size(), DisplacementDim);
        }

        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    ThermoMechanicalPhaseFieldProcessData<DisplacementDim> process_data{
        materialIDs(mesh),
        std::move(solid_constitutive_relations),
        residual_stiffness,
        crack_resistance,
        crack_length_scale,
        kinetic_coefficient,
        solid_density,
        linear_thermal_expansion_coefficient,
        specific_heat_capacity,
        thermal_conductivity,
        residual_thermal_conductivity,
        specific_body_force,
        reference_temperature};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<ThermoMechanicalPhaseFieldProcess<DisplacementDim>>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        mechanics_related_process_id, phase_field_process_id,
        heat_conduction_process_id);
}

template std::unique_ptr<Process> createThermoMechanicalPhaseFieldProcess<2>(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<Process> createThermoMechanicalPhaseFieldProcess<3>(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace ThermoMechanicalPhaseField
}  // namespace ProcessLib
