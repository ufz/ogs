/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateThermoMechanicalPhaseFieldProcess.h"

#include <cassert>

#include "MaterialLib/SolidModels/CreateLinearElasticIsotropic.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropicPhaseField.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"

#include "ThermoMechanicalPhaseFieldProcess.h"
#include "ThermoMechanicalPhaseFieldProcessData.h"

namespace ProcessLib
{
namespace ThermoMechanicalPhaseField
{
template <int DisplacementDim>
std::unique_ptr<Process> createThermoMechanicalPhaseFieldProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "THERMO_MECHANICAL_PHASE_FIELD");
    DBUG("Create ThermoMechanicalPhaseFieldProcess.");

    auto const staggered_scheme =
        //! \ogs_file_param{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(staggered_scheme && (*staggered_scheme == "staggered"));
    if (!(staggered_scheme && (*staggered_scheme == "staggered")))
    {
        OGS_FATAL(
            "Monolithic scheme for thermo-mechanical coupled "
            "phase field process is currently not available.");
    }

    // Process variable.

    //! \ogs_file_param{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");
    ProcessVariable* variable_T;
    ProcessVariable* variable_ph;
    ProcessVariable* variable_u;
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    int mechanics_related_process_id = 0;
    int phase_field_process_id = 0;
    int heat_conduction_process_id = 0;
    if (use_monolithic_scheme)  // monolithic scheme.
    {
        auto per_process_variables = findProcessVariables(
            variables, pv_config,
            {//! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__process_variables__temperature}
             "temperature",
             //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__process_variables__phasefield}
             "phasefield",
             //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__process_variables__displacement}
             "displacement"});
        variable_T = &per_process_variables[0].get();
        variable_ph = &per_process_variables[1].get();
        variable_u = &per_process_variables[2].get();
        process_variables.push_back(std::move(per_process_variables));
    }
    else  // staggered scheme.
    {
        using namespace std::string_literals;
        for (auto const& variable_name :
             {"temperature"s, "displacement"s, "phasefield"s})
        {
            auto per_process_variables =
                findProcessVariables(variables, pv_config, {variable_name});
            process_variables.push_back(std::move(per_process_variables));
        }
        heat_conduction_process_id = 0;
        mechanics_related_process_id = 1;
        phase_field_process_id = 2;

        variable_T = &process_variables[heat_conduction_process_id][0].get();
        variable_u = &process_variables[mechanics_related_process_id][0].get();
        variable_ph = &process_variables[phase_field_process_id][0].get();
    }

    DBUG("Associate displacement with process variable \'%s\'.",
         variable_u->getName().c_str());

    if (variable_u->getNumberOfComponents() != DisplacementDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '%s' is different "
            "from the displacement dimension: got %d, expected %d",
            variable_u->getName().c_str(),
            variable_u->getNumberOfComponents(),
            DisplacementDim);
    }

    DBUG("Associate phase field with process variable \'%s\'.",
         variable_ph->getName().c_str());
    if (variable_ph->getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Phase field process variable '%s' is not a scalar variable but "
            "has "
            "%d components.",
            variable_ph->getName().c_str(),
            variable_ph->getNumberOfComponents());
    }

    DBUG("Associate temperature with process variable \'%s\'.",
         variable_T->getName().c_str());
    if (variable_T->getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Temperature process variable '%s' is not a scalar variable but "
            "has "
            "%d components.",
            variable_T->getName().c_str(),
            variable_T->getNumberOfComponents());
    }

    // Constitutive relation.
    // read type;
    auto const constitutive_relation_config =
        //! \ogs_file_param{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__constitutive_relation}
        config.getConfigSubtree("constitutive_relation");

    auto const type =
        //! \ogs_file_param{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__constitutive_relation__type}
        constitutive_relation_config.peekConfigParameter<std::string>("type");

    std::unique_ptr<MaterialLib::Solids::PhaseFieldExtension<DisplacementDim>>
        material = nullptr;
    if (type == "LinearElasticIsotropic")
    {
        auto elastic_model =
            MaterialLib::Solids::createLinearElasticIsotropic<DisplacementDim>(
                parameters, constitutive_relation_config);
        material = std::make_unique<
            MaterialLib::Solids::LinearElasticIsotropicPhaseField<
                DisplacementDim>>(
            std::move(elastic_model->getMaterialProperties()));
    }
    else
    {
        OGS_FATAL(
            "Cannot construct constitutive relation of given type \'%s\'.",
            type.c_str());
    }

    // Residual stiffness
    auto& residual_stiffness = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__residual_stiffness}
        "residual_stiffness", parameters, 1);
    DBUG("Use \'%s\' as residual stiffness.", residual_stiffness.name.c_str());

    // Crack resistance
    auto& crack_resistance = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__crack_resistance}
        "crack_resistance", parameters, 1);
    DBUG("Use \'%s\' as crack resistance.", crack_resistance.name.c_str());

    // Crack length scale
    auto& crack_length_scale = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__crack_length_scale}
        "crack_length_scale", parameters, 1);
    DBUG("Use \'%s\' as crack length scale.", crack_length_scale.name.c_str());

    // Kinetic coefficient
    auto& kinetic_coefficient = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__kinetic_coefficient}
        "kinetic_coefficient", parameters, 1);
    DBUG("Use \'%s\' as kinetic coefficient.",
         kinetic_coefficient.name.c_str());

    // Solid density
    auto& solid_density = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__solid_density}
        "solid_density", parameters, 1);
    DBUG("Use \'%s\' as solid density parameter.", solid_density.name.c_str());

    // History field
    auto& history_field = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__history_field}
        "history_field", parameters, 1);
    DBUG("Use \'%s\' as history field.", history_field.name.c_str());

    // Linear thermal expansion coefficient
    auto& linear_thermal_expansion_coefficient = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__linear_thermal_expansion_coefficient}
        "linear_thermal_expansion_coefficient", parameters, 1);
    DBUG("Use \'%s\' as linear thermal expansion coefficient.",
         linear_thermal_expansion_coefficient.name.c_str());

    // Specific heat capacity
    auto& specific_heat_capacity = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__specific_heat_capacity}
        "specific_heat_capacity", parameters, 1);
    DBUG("Use \'%s\' as specific heat capacity.",
         specific_heat_capacity.name.c_str());

    // Thermal conductivity
    auto& thermal_conductivity = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__thermal_conductivity}
        "thermal_conductivity", parameters, 1);
    DBUG("Use \'%s\' as thermal conductivity parameter.",
         thermal_conductivity.name.c_str());
    // Residual thermal conductivity
    auto& residual_thermal_conductivity = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__residual_thermal_conductivity}
        "residual_thermal_conductivity", parameters, 1);
    DBUG("Use \'%s\' as residual thermal conductivity parameter.",
         residual_thermal_conductivity.name.c_str());
    // Reference temperature
    const double reference_temperature =
        //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__reference_temperature}
        config.getConfigParameter<double>("reference_temperature");

    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__THERMO_MECHANICAL_PHASE_FIELD__specific_body_force}
            config.getConfigParameter<std::vector<double>>(
                "specific_body_force");
        if (specific_body_force.size() != DisplacementDim)
            OGS_FATAL(
                "The size of the specific body force vector does not match the "
                "displacement dimension. Vector size is %d, displacement "
                "dimension is %d",
                specific_body_force.size(), DisplacementDim);

        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    ThermoMechanicalPhaseFieldProcessData<DisplacementDim> process_data{
        std::move(material),
        residual_stiffness,
        crack_resistance,
        crack_length_scale,
        kinetic_coefficient,
        solid_density,
        history_field,
        linear_thermal_expansion_coefficient,
        specific_heat_capacity,
        thermal_conductivity,
        residual_thermal_conductivity,
        reference_temperature,
        specific_body_force};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"temperature_phasefield_displacement"});

    ProcessLib::createSecondaryVariables(config, secondary_variables,
                                         named_function_caller);

    return std::make_unique<ThermoMechanicalPhaseFieldProcess<DisplacementDim>>(
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller),
        use_monolithic_scheme, mechanics_related_process_id,
        phase_field_process_id, heat_conduction_process_id);
}

template std::unique_ptr<Process> createThermoMechanicalPhaseFieldProcess<2>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<Process> createThermoMechanicalPhaseFieldProcess<3>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace ThermoMechanicalPhaseField
}  // namespace ProcessLib
