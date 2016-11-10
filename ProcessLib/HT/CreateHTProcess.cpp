/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateHTProcess.h"

#include "HTProcess.h"
#include "HTProcessData.h"
#include "ProcessLib/Parameter/ConstantParameter.h"
#include "ProcessLib/Utils/ParseSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
namespace HT
{
std::unique_ptr<Process> createHTProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{process__type}
    config.checkConfigParameter("type", "HT");

    DBUG("Create HTProcess.");

    // Process variable.
    auto process_variables = findProcessVariables(
        variables, config,
        {"temperature", "pressure"});  // configure two Pcs

    // Porosity parameter.
    auto& porosity = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__HT__porosity}
        "porosity", parameters, 1);
    DBUG("Use \'%s\' as porosity parameter.", porosity.name.c_str());

    // Parameter for the intrinsic permeability (only one scalar per element,
    // i.e., the isotropic case is handled at the moment)
    auto& intrinsic_permeability = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__HT__intrinsic_permeability}
        "intrinsic_permeability", parameters, 1);
    DBUG("Use \'%s\' as intrinsic_permeability parameter.",
         intrinsic_permeability.name.c_str());

    // Parameter for the specific storage.
    auto& specific_storage = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__HT__specific_storage}
        "specific_storage", parameters, 1);
    DBUG("Use \'%s\' as specific storage parameter.", specific_storage.name.c_str());

    // Parameter for the reference_temperature.
    auto& reference_temperature_fluid_density_model = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__HT__reference_temperature}
        "reference_temperature_fluid_density_model", parameters, 1);
    DBUG(
        "Use \'%s\' as reference temperature for the fluid density model "
        "parameter.",
        reference_temperature_fluid_density_model.name.c_str());

    // Parameter for the viscosity.
    auto& viscosity = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__HT__viscosity}
        "viscosity", parameters, 1);
    DBUG("Use \'%s\' as viscosity parameter.", viscosity.name.c_str());

    // Parameter for the density of the solid.
    auto& density_solid = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__HT__density_solid}
        "density_solid", parameters, 1);
    DBUG("Use \'%s\' as density_solid parameter.", density_solid.name.c_str());

    // Parameter for the density of the fluid.
    auto& density_fluid = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__HT__density_fluid}
        "density_fluid", parameters, 1);
    DBUG("Use \'%s\' as density_fluid parameter.", density_fluid.name.c_str());

    // Parameter for the specific heat capacity of the solid.
    auto& specific_heat_capacity_solid = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__HT__specific_heat_capacity_solid}
        "specific_heat_capacity_solid", parameters, 1);
    DBUG("Use \'%s\' as specific_heat_capacity_solid parameter.",
         specific_heat_capacity_solid.name.c_str());

    // Parameter for the specific heat capacity of the fluid.
    auto& specific_heat_capacity_fluid = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__HT__specific_heat_capacity_fluid}
        "specific_heat_capacity_fluid", parameters, 1);
    DBUG("Use \'%s\' as specific_heat_capacity_fluid parameter.",
         specific_heat_capacity_fluid.name.c_str());

    // Parameter for the thermal conductivity of the solid (only one scalar per
    // element, i.e., the isotropic case is handled at the moment)
    auto& thermal_conductivity_solid = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__HT__thermal_conductivity_solid}
        "thermal_conductivity_solid", parameters, 1);
    DBUG("Use \'%s\' as thermal_conductivity_solid parameter.",
         thermal_conductivity_solid.name.c_str());

    // Parameter for the thermal conductivity of the fluid.
    auto& thermal_conductivity_fluid = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__HT__thermal_conductivity_fluid}
        "thermal_conductivity_fluid", parameters, 1);
    DBUG("Use \'%s\' as thermal_conductivity_fluid parameter.",
         thermal_conductivity_fluid.name.c_str());

    // Parameter for the thermal expansion coefficient.
    auto& thermal_expansion_coefficient_fluid = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__HT__thermal_expansion_coefficient_fluid}
        "thermal_expansion_coefficient_fluid", parameters, 1);
    DBUG("Use \'%s\' as thermal expansion coefficient of the fluid parameter.",
         thermal_expansion_coefficient_fluid.name.c_str());

    // Specific body force parameter.
    Eigen::Vector3d specific_body_force;
    std::vector<double> const b =
        //! \ogs_file_param_special{process__HT__specific_body_force}
        config.getConfigParameter<std::vector<double>>("specific_body_force");
    assert(b.size() > 0 && b.size() < 4);
    bool const has_gravity = MathLib::toVector(b).norm() > 0;
    if (has_gravity)
        std::copy_n(b.data(), b.size(), specific_body_force.data());

    HTProcessData process_data{
        porosity,
        intrinsic_permeability,
        specific_storage,
        viscosity,
        density_solid,
        density_fluid,
        specific_heat_capacity_solid,
        specific_heat_capacity_fluid,
        thermal_conductivity_solid,
        thermal_conductivity_fluid,
        thermal_expansion_coefficient_fluid,
        reference_temperature_fluid_density_model,
        specific_body_force,
        has_gravity};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"HT_temperature_pressure"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    return std::unique_ptr<Process>{new HTProcess{
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller)}};
}

}  // namespace HT
}  // namespace ProcessLib
