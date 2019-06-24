/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateThermoHydroMechanicsProcess.h"

#include <cassert>

#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "ThermoHydroMechanicsProcess.h"
#include "ThermoHydroMechanicsProcessData.h"

namespace ProcessLib
{
namespace ThermoHydroMechanics
{
template <int DisplacementDim>
std::unique_ptr<Process> createThermoHydroMechanicsProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "THERMO_HYDRO_MECHANICS");
    DBUG("Create ThermoHydroMechanicsProcess.");

    auto const staggered_scheme =
        //! \ogs_file_param{prj__processes__process__THERMO_HYDRO_MECHANICS__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(staggered_scheme && (*staggered_scheme == "staggered"));

    // Process variable.

    //! \ogs_file_param{prj__processes__process__THERMO_HYDRO_MECHANICS__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    ProcessVariable* variable_T;
    ProcessVariable* variable_p;
    ProcessVariable* variable_u;
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    if (use_monolithic_scheme)  // monolithic scheme.
    {
        auto per_process_variables = findProcessVariables(
            variables, pv_config,
            {//! \ogs_file_param_special{prj__processes__process__THERMO_HYDRO_MECHANICS__process_variables__temperature}
             "temperature",
             //! \ogs_file_param_special{prj__processes__process__THERMO_HYDRO_MECHANICS__process_variables__pressure}
             "pressure",
             //! \ogs_file_param_special{prj__processes__process__THERMO_HYDRO_MECHANICS__process_variables__displacement}
             "displacement"});
        variable_T = &per_process_variables[0].get();
        variable_p = &per_process_variables[1].get();
        variable_u = &per_process_variables[2].get();
        process_variables.push_back(std::move(per_process_variables));
    }
    else  // staggered scheme.
    {
        using namespace std::string_literals;
        for (auto const& variable_name :
             {"temperature"s, "pressure"s, "displacement"s})
        {
            auto per_process_variables =
                findProcessVariables(variables, pv_config, {variable_name});
            process_variables.push_back(std::move(per_process_variables));
        }
        variable_T = &process_variables[0][0].get();
        variable_p = &process_variables[1][0].get();
        variable_u = &process_variables[2][0].get();
    }

    DBUG("Associate displacement with process variable '%s'.",
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

    DBUG("Associate pressure with process variable '%s'.",
         variable_p->getName().c_str());
    if (variable_p->getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Pressure process variable '%s' is not a scalar variable but has "
            "%d components.",
            variable_p->getName().c_str(),
            variable_p->getNumberOfComponents());
    }

    DBUG("Associate temperature with process variable '%s'.",
         variable_T->getName().c_str());
    if (variable_T->getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "temperature process variable '%s' is not a scalar variable but "
            "has %d components.",
            variable_T->getName().c_str(),
            variable_T->getNumberOfComponents());
    }

    auto solid_constitutive_relations =
        MaterialLib::Solids::createConstitutiveRelations<DisplacementDim>(
            parameters, local_coordinate_system, config);

    // Intrinsic permeability (only one scalar per element, i.e. the isotropic
    // case is handled at the moment)
    auto& intrinsic_permeability = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_HYDRO_MECHANICS__intrinsic_permeability}
        "intrinsic_permeability", parameters, 1, &mesh);

    DBUG("Use '%s' as intrinsic conductivity parameter.",
         intrinsic_permeability.name.c_str());

    // Storage coefficient
    auto& specific_storage = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_HYDRO_MECHANICS__specific_storage}
        "specific_storage", parameters, 1, &mesh);

    DBUG("Use '%s' as storage coefficient parameter.",
         specific_storage.name.c_str());

    // Fluid viscosity
    auto& fluid_viscosity = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_HYDRO_MECHANICS__fluid_viscosity}
        "fluid_viscosity", parameters, 1, &mesh);
    DBUG("Use '%s' as fluid viscosity parameter.",
         fluid_viscosity.name.c_str());

    // Fluid density
    auto& fluid_density = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_HYDRO_MECHANICS__fluid_density}
        "fluid_density", parameters, 1, &mesh);
    DBUG("Use '%s' as fluid density parameter.", fluid_density.name.c_str());

    // Biot coefficient
    auto& biot_coefficient = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_HYDRO_MECHANICS__biot_coefficient}
        "biot_coefficient", parameters, 1, &mesh);
    DBUG("Use '%s' as Biot coefficient parameter.",
         biot_coefficient.name.c_str());

    // Porosity
    auto& porosity = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_HYDRO_MECHANICS__porosity}
        "porosity", parameters, 1, &mesh);
    DBUG("Use '%s' as porosity parameter.", porosity.name.c_str());

    // Solid density
    auto& solid_density = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_HYDRO_MECHANICS__solid_density}
        "solid_density", parameters, 1, &mesh);
    DBUG("Use '%s' as solid density parameter.", solid_density.name.c_str());

    // linear thermal expansion coefficient for solid
    auto const& solid_linear_thermal_expansion_coefficient =
        ParameterLib::findParameter<double>(
            config,
            //! \ogs_file_param_special{prj__processes__process__THERMO_HYDRO_MECHANICS__solid_linear_thermal_expansion_coefficient}
            "solid_linear_thermal_expansion_coefficient", parameters, 1, &mesh);
    DBUG("Use '%s' as solid linear thermal expansion coefficient parameter.",
         solid_linear_thermal_expansion_coefficient.name.c_str());

    // volumetric thermal expansion coefficient for fluid
    auto const& fluid_volumetric_thermal_expansion_coefficient =
        ParameterLib::findParameter<double>(
            config,
            //! \ogs_file_param_special{prj__processes__process__THERMO_HYDRO_MECHANICS__fluid_volumetric_thermal_expansion_coefficient}
            "fluid_volumetric_thermal_expansion_coefficient", parameters, 1,
            &mesh);
    DBUG(
        "Use '%s' as fluid volumetric thermal expansion coefficient "
        "parameter.",
        fluid_volumetric_thermal_expansion_coefficient.name.c_str());

    // specific heat capacity for solid
    auto& solid_specific_heat_capacity = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_HYDRO_MECHANICS__solid_specific_heat_capacity}
        "solid_specific_heat_capacity", parameters, 1, &mesh);
    DBUG("Use '%s' as solid specific heat capacity parameter.",
         solid_specific_heat_capacity.name.c_str());

    // specific heat capacity for fluid
    auto& fluid_specific_heat_capacity = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_HYDRO_MECHANICS__fluid_specific_heat_capacity}
        "fluid_specific_heat_capacity", parameters, 1, &mesh);
    DBUG("Use '%s' as fluid specific heat capacity parameter.",
         fluid_specific_heat_capacity.name.c_str());

    // thermal conductivity for solid // currently only considers isotropic
    auto& solid_thermal_conductivity = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_HYDRO_MECHANICS__solid_thermal_conductivity}
        "solid_thermal_conductivity", parameters, 1, &mesh);
    DBUG("Use '%s' as solid thermal conductivity parameter.",
         solid_thermal_conductivity.name.c_str());

    // thermal conductivity for fluid // currently only considers isotropic
    auto& fluid_thermal_conductivity = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_HYDRO_MECHANICS__fluid_thermal_conductivity}
        "fluid_thermal_conductivity", parameters, 1, &mesh);
    DBUG("Use '%s' as fluid thermal conductivity parameter.",
         fluid_thermal_conductivity.name.c_str());

    // reference temperature
    auto& reference_temperature = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_HYDRO_MECHANICS__reference_temperature}
        "reference_temperature", parameters, 1, &mesh);
    DBUG("Use '%s' as reference temperature parameter.",
         reference_temperature.name.c_str());

    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__THERMO_HYDRO_MECHANICS__specific_body_force}
            config.getConfigParameter<std::vector<double>>(
                "specific_body_force");
        if (b.size() != DisplacementDim)
        {
            OGS_FATAL(
                "The size of the specific body force vector does not match the "
                "displacement dimension. Vector size is %d, displacement "
                "dimension is %d",
                b.size(), DisplacementDim);
        }

        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    ThermoHydroMechanicsProcessData<DisplacementDim> process_data{
        materialIDs(mesh),
        std::move(solid_constitutive_relations),
        intrinsic_permeability,
        specific_storage,
        fluid_viscosity,
        fluid_density,
        biot_coefficient,
        porosity,
        solid_density,
        solid_linear_thermal_expansion_coefficient,
        fluid_volumetric_thermal_expansion_coefficient,
        solid_specific_heat_capacity,
        fluid_specific_heat_capacity,
        solid_thermal_conductivity,
        fluid_thermal_conductivity,
        reference_temperature,
        specific_body_force};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"ThermoHydroMechanics_displacement"});

    ProcessLib::createSecondaryVariables(config, secondary_variables,
                                         named_function_caller);

    return std::make_unique<ThermoHydroMechanicsProcess<DisplacementDim>>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        std::move(named_function_caller), use_monolithic_scheme);
}

template std::unique_ptr<Process> createThermoHydroMechanicsProcess<2>(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<Process> createThermoHydroMechanicsProcess<3>(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace ThermoHydroMechanics
}  // namespace ProcessLib
