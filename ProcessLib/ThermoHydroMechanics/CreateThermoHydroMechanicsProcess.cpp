/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateThermoHydroMechanicsProcess.h"

#include <cassert>

#include "MaterialLib/SolidModels/CreateLinearElasticIsotropic.h"
#include "ProcessLib/Utils/ParseSecondaryVariables.h"

#include "ThermoHydroMechanicsProcess.h"
#include "ThermoHydroMechanicsProcessData.h"

namespace ProcessLib
{
namespace ThermoHydroMechanics
{
template <int DisplacementDim>
class ThermoHydroMechanicsProcess;

extern template class ThermoHydroMechanicsProcess<2>;
extern template class ThermoHydroMechanicsProcess<3>;

template <int DisplacementDim>
std::unique_ptr<Process> createThermoHydroMechanicsProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{process__type}
    config.checkConfigParameter("type", "THERMO_HYDRO_MECHANICS");
    DBUG("Create ThermoHydroMechanicsProcess.");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    auto process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{process__THERMO_HYDRO_MECHANICS_process_variables__process_variable}
          "temperature","pressure", "displacement"});

    DBUG("Associate displacement with process variable \'%s\'.",
         process_variables[2].get().getName().c_str());

    if (process_variables[2].get().getNumberOfComponents() !=
        DisplacementDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '%s' is different "
            "from the displacement dimension: got %d, expected %d",
            process_variables[2].get().getName().c_str(),
            process_variables[2].get().getNumberOfComponents(),
            DisplacementDim);
    }

    DBUG("Associate pressure with process variable \'%s\'.",
         process_variables[1].get().getName().c_str());
    if (process_variables[1].get().getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Pressure process variable '%s' is not a scalar variable but has "
            "%d components.",
            process_variables[1].get().getName().c_str(),
            process_variables[1].get().getNumberOfComponents(),
            DisplacementDim);
    }

    DBUG("Associate temperature with process variable \'%s\'.",
         process_variables[0].get().getName().c_str());
    if (process_variables[0].get().getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Pressure process variable '%s' is not a scalar variable but has "
            "%d components.",
            process_variables[0].get().getName().c_str(),
            process_variables[0].get().getNumberOfComponents(),
            DisplacementDim);
    }

    // Constitutive relation.
    // read type;
    auto const constitutive_relation_config =
        //! \ogs_file_param{process__THERMO_HYDRO_MECHANICS_constitutive_relation}
        config.getConfigSubtree("constitutive_relation");

    auto const type =
        constitutive_relation_config.peekConfigParameter<std::string>("type");

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        material = nullptr;
    if (type == "LinearElasticIsotropic")
    {
        material =
            MaterialLib::Solids::createLinearElasticIsotropic<DisplacementDim>(
                parameters, constitutive_relation_config);
    }
    else
    {
        OGS_FATAL(
            "Cannot construct constitutive relation of given type \'%s\'.",
            type.c_str());
    }

    // Intrinsic permeability
    auto& intrinsic_permeability = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__THERMO_HYDRO_MECHANICS_intrinsic_permeability}
        "intrinsic_permeability",
        parameters, 1);

    DBUG("Use \'%s\' as hydraulic conductivity parameter.",
         intrinsic_permeability.name.c_str());

    // Storage coefficient
    auto& storage_coefficient = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__THERMO_HYDRO_MECHANICS_storage_coefficient}
        "storage_coefficient",
        parameters, 1);

    DBUG("Use \'%s\' as storage coefficient parameter.",
         storage_coefficient.name.c_str());


    // Fluid viscosity
    auto& fluid_viscosity = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__THERMO_HYDRO_MECHANICS_fluid_viscosity}
        "fluid_viscosity",
        parameters, 1);
    DBUG("Use \'%s\' as fluid viscosity parameter.",
         fluid_viscosity.name.c_str());

    // Biot coefficient
    auto& biot_coefficient = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__THERMO_HYDRO_MECHANICS_biot_coefficient}
        "biot_coefficient",
        parameters, 1);
    DBUG("Use \'%s\' as Biot coefficient parameter.",
         biot_coefficient.name.c_str());

    // Porosity
    auto& porosity = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__THERMO_HYDRO_MECHANICS_porosity}
        "porosity",
        parameters, 1);
    DBUG("Use \'%s\' as porosity parameter.",
         porosity.name.c_str());

    // Solid density
    auto& solid_density = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__THERMO_HYDRO_MECHANICS_solid_density}
        "solid_density",
        parameters, 1);
    DBUG("Use \'%s\' as solid density parameter.",
         solid_density.name.c_str());

    // Fluid density
    auto const& fluid_density = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__THERMO_HYDRO_MECHANICS_fluid_density}
        "fluid_density",
        parameters, 1);
    DBUG("Use \'%s\' as fluid density parameter.",
         fluid_density.name.c_str());

    // thermal expansion coefficient for solid
    auto const& beta_solid = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__THERMO_HYDRO_MECHANICS_solid_thermal_expansion_coefficient}
        "beta_solid",
        parameters, 1);
    DBUG("Use \'%s\' as solid thermal expansion coefficient parameter.",
         beta_solid.name.c_str());

    // thermal expansion coefficient for fluid
    auto& beta_fluid = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__THERMO_HYDRO_MECHANICS_fluid_thermal_expansion_coefficient}
        "beta_fluid",
        parameters, 1);
    DBUG("Use \'%s\' as fluid thermal expansion coefficient parameter.",
         beta_fluid.name.c_str());

    // specific heat capacity for fluid
    auto& fluid_heat_capacity = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__THERMO_HYDRO_MECHANICS_fluid_heat_capacity}
        "fluid_heat_capacity",
        parameters, 1);
    DBUG("Use \'%s\' as fluid heat capacity parameter.",
         fluid_heat_capacity.name.c_str());

    // specific heat capacity for solid
    auto& solid_heat_capacity = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__THERMO_HYDRO_MECHANICS_solid_heat_capacity}
        "solid_heat_capacity",
        parameters, 1);
    DBUG("Use \'%s\' as solid heat capacity parameter.",
         solid_heat_capacity.name.c_str());

    // thermal conductivity for solid
    auto& lambda_s = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__THERMO_HYDRO_MECHANICS_solid_thermal_conductivity}
        "lambda_s",
        parameters, 1);
    DBUG("Use \'%s\' as solid thermal conductivity parameter.",
         lambda_s.name.c_str());

    // thermal conductivity for fluid
    auto& lambda_f = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__THERMO_HYDRO_MECHANICS_fluid_thermal_conductivity}
        "lambda_f",
        parameters, 1);
    DBUG("Use \'%s\' as fluid thermal conductivity parameter.",
         lambda_f.name.c_str());

    // reference temperature
    auto& reference_temperature = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__THERMO_HYDRO_MECHANICS_reference_temperature}
        "reference_temperature",
        parameters, 1);
    DBUG("Use \'%s\' as reference temperature parameter.",
         reference_temperature.name.c_str());

    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__specific_body_force}
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

    ThermoHydroMechanicsProcessData<DisplacementDim> process_data{
        std::move(material), intrinsic_permeability, storage_coefficient,
        fluid_viscosity,     biot_coefficient,       porosity,
        solid_density,       fluid_density,          beta_solid,
        beta_fluid,          fluid_heat_capacity,    solid_heat_capacity,
        lambda_s,            lambda_f,            reference_temperature,
        specific_body_force
    };

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"ThermoHydroMechanics_displacement"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    return std::unique_ptr<ThermoHydroMechanicsProcess<DisplacementDim>>{
        new ThermoHydroMechanicsProcess<DisplacementDim>{
            mesh, std::move(jacobian_assembler), parameters, integration_order,
            std::move(process_variables), std::move(process_data),
            std::move(secondary_variables), std::move(named_function_caller)}};
}

template std::unique_ptr<Process> createThermoHydroMechanicsProcess<2>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<Process> createThermoHydroMechanicsProcess<3>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace ThermoHydroMechanics
}  // namespace ProcessLib
