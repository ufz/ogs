/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateThermoMechanicsProcess.h"

#include <cassert>

#include "MaterialLib/SolidModels/CreateEhlers.h"
#include "MaterialLib/SolidModels/CreateLinearElasticIsotropic.h"
#include "MaterialLib/SolidModels/CreateLubby2.h"
#include "ProcessLib/Utils/ParseSecondaryVariables.h"

#include "ThermoMechanicsProcess.h"
#include "ThermoMechanicsProcessData.h"

namespace ProcessLib
{
namespace ThermoMechanics
{
template <int DisplacementDim>
class ThermoMechanicsProcess;

extern template class ThermoMechanicsProcess<2>;
extern template class ThermoMechanicsProcess<3>;

template <int DisplacementDim>
std::unique_ptr<Process> createThermoMechanicsProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "THERMO_MECHANICS");
    DBUG("Create ThermoMechanicsProcess.");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__THERMO_MECHANICS__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    auto process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICS__process_variables__temperature}
         "temperature",
         //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICS__process_variables__displacement}
         "displacement"});

    DBUG("Associate displacement with process variable \'%s\'.",
         process_variables[1].get().getName().c_str());

    if (process_variables[1].get().getNumberOfComponents() != DisplacementDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '%s' is different "
            "from the displacement dimension: got %d, expected %d",
            process_variables[1].get().getName().c_str(),
            process_variables[1].get().getNumberOfComponents(),
            DisplacementDim);
    }

    DBUG("Associate temperature with process variable \'%s\'.",
         process_variables[0].get().getName().c_str());
    if (process_variables[0].get().getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Temperature process variable '%s' is not a scalar variable but has "
            "%d components.",
            process_variables[0].get().getName().c_str(),
            process_variables[0].get().getNumberOfComponents(),
            DisplacementDim);
    }

    // Constitutive relation.
    // read type;
    auto const constitutive_relation_config =
        //! \ogs_file_param{prj__processes__process__THERMO_MECHANICS__constitutive_relation}
        config.getConfigSubtree("constitutive_relation");

    auto const type =
        //! \ogs_file_param{prj__processes__process__THERMO_MECHANICS__constitutive_relation__type}
        constitutive_relation_config.peekConfigParameter<std::string>("type");

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        material = nullptr;
    if (type == "Ehlers")
    {
        material = MaterialLib::Solids::Ehlers::createEhlers<DisplacementDim>(
            parameters, constitutive_relation_config);
    }
    else if (type == "LinearElasticIsotropic")
    {
        material =
            MaterialLib::Solids::createLinearElasticIsotropic<DisplacementDim>(
                parameters, constitutive_relation_config);
    }
    else if (type == "Lubby2")
    {
        material = MaterialLib::Solids::Lubby2::createLubby2<DisplacementDim>(
            parameters, constitutive_relation_config);
    }
    else
    {
        OGS_FATAL(
            "Cannot construct constitutive relation of given type \'%s\'.",
            type.c_str());
    }

    // Reference solid density
    auto& reference_solid_density = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICS__reference_solid_density}
        "reference_solid_density", parameters, 1);
    DBUG("Use \'%s\' as solid density parameter.",
         reference_solid_density.name.c_str());

    // Linear thermal expansion coefficient
    auto& linear_thermal_expansion_coefficient = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICS__linear_thermal_expansion_coefficient}
        "linear_thermal_expansion_coefficient", parameters, 1);
    DBUG("Use \'%s\' as linear thermal expansion coefficient.",
         linear_thermal_expansion_coefficient.name.c_str());
    // Specific heat capacity
    auto& specific_heat_capacity = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICS__specific_heat_capacity}
        "specific_heat_capacity", parameters, 1);
    DBUG("Use \'%s\' as specific heat capacity parameter.",
         specific_heat_capacity.name.c_str());
    // Thermal conductivity // TODO To be changed as tensor input.
    auto& thermal_conductivity = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_MECHANICS__thermal_conductivity}
        "thermal_conductivity", parameters, 1);
    DBUG("Use \'%s\' as thermal conductivity parameter.",
         thermal_conductivity.name.c_str());
    // Reference temperature
    const double reference_temperature =
        //! \ogs_file_param{prj__processes__process__THERMO_MECHANICS__reference_temperature}
        config.getConfigParameter<double>("reference_temperature");

    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__THERMO_MECHANICS__specific_body_force}
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

    ThermoMechanicsProcessData<DisplacementDim> process_data{
        std::move(material),
        reference_solid_density,
        linear_thermal_expansion_coefficient,
        specific_heat_capacity,
        thermal_conductivity,
        reference_temperature,
        specific_body_force};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"ThermoMechanics_temperature_displacement"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    return std::unique_ptr<ThermoMechanicsProcess<DisplacementDim>>{
        new ThermoMechanicsProcess<DisplacementDim>{
            mesh, std::move(jacobian_assembler), parameters, integration_order,
            std::move(process_variables), std::move(process_data),
            std::move(secondary_variables), std::move(named_function_caller)}};
}

template std::unique_ptr<Process> createThermoMechanicsProcess<2>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<Process> createThermoMechanicsProcess<3>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace ThermoMechanics
}  // namespace ProcessLib
