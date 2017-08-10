/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateHydroMechanicsProcess.h"

#include <cassert>

#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MaterialLib/SolidModels/CreateLinearElasticIsotropic.h"
#include "ProcessLib/Utils/ParseSecondaryVariables.h"

#include "HydroMechanicsProcess.h"
#include "HydroMechanicsProcessData.h"

namespace ProcessLib
{
namespace HydroMechanics
{

template <int DisplacementDim>
std::unique_ptr<Process> createHydroMechanicsProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "HYDRO_MECHANICS");
    DBUG("Create HydroMechanicsProcess.");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    auto process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS__process_variables__pressure}
         "pressure",
         //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS__process_variables__displacement}
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

    DBUG("Associate pressure with process variable \'%s\'.",
         process_variables[0].get().getName().c_str());
    if (process_variables[0].get().getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Pressure process variable '%s' is not a scalar variable but has "
            "%d components.",
            process_variables[0].get().getName().c_str(),
            process_variables[0].get().getNumberOfComponents());
    }

    // Constitutive relation.
    // read type;
    auto const constitutive_relation_config =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__constitutive_relation}
        config.getConfigSubtree("constitutive_relation");

    auto const type =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__constitutive_relation__type}
        constitutive_relation_config.peekConfigParameter<std::string>("type");

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        material = nullptr;
    if (type == "LinearElasticIsotropic")
    {
        material =
            MaterialLib::Solids::createLinearElasticIsotropic<DisplacementDim,
                MaterialLib::Solids::LinearElasticIsotropic>(
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
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS__intrinsic_permeability}
        "intrinsic_permeability", parameters, 1);

    DBUG("Use \'%s\' as intrinsic conductivity parameter.",
         intrinsic_permeability.name.c_str());

    // Storage coefficient
    auto& specific_storage = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS__specific_storage}
        "specific_storage", parameters, 1);

    DBUG("Use \'%s\' as storage coefficient parameter.",
         specific_storage.name.c_str());

    // Fluid viscosity
    auto& fluid_viscosity = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS__fluid_viscosity}
        "fluid_viscosity", parameters, 1);
    DBUG("Use \'%s\' as fluid viscosity parameter.",
         fluid_viscosity.name.c_str());

    // Fluid density
    auto& fluid_density = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS__fluid_density}
        "fluid_density", parameters, 1);
    DBUG("Use \'%s\' as fluid density parameter.", fluid_density.name.c_str());

    // Biot coefficient
    auto& biot_coefficient = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS__biot_coefficient}
        "biot_coefficient", parameters, 1);
    DBUG("Use \'%s\' as Biot coefficient parameter.",
         biot_coefficient.name.c_str());

    // Porosity
    auto& porosity = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS__porosity}
        "porosity", parameters, 1);
    DBUG("Use \'%s\' as porosity parameter.", porosity.name.c_str());

    // Solid density
    auto& solid_density = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HYDRO_MECHANICS__solid_density}
        "solid_density", parameters, 1);
    DBUG("Use \'%s\' as solid density parameter.", solid_density.name.c_str());

    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__specific_body_force}
            config.getConfigParameter<std::vector<double>>(
                "specific_body_force");
        if (b.size() != DisplacementDim)
            OGS_FATAL(
                "The size of the specific body force vector does not match the "
                "displacement dimension. Vector size is %d, displacement "
                "dimension is %d",
                b.size(), DisplacementDim);

        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    HydroMechanicsProcessData<DisplacementDim> process_data{
        std::move(material),
        intrinsic_permeability,
        specific_storage,
        fluid_viscosity,
        fluid_density,
        biot_coefficient,
        porosity,
        solid_density,
        specific_body_force};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"HydroMechanics_displacement"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    return std::make_unique<HydroMechanicsProcess<DisplacementDim>>(
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller));
}

template std::unique_ptr<Process> createHydroMechanicsProcess<2>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<Process> createHydroMechanicsProcess<3>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace HydroMechanics
}  // namespace ProcessLib
