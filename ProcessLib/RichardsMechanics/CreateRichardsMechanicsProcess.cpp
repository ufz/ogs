/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateRichardsMechanicsProcess.h"

#include <cassert>

#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/RichardsFlow/CreateRichardsFlowMaterialProperties.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "RichardsMechanicsProcess.h"
#include "RichardsMechanicsProcessData.h"

namespace ProcessLib
{
namespace RichardsMechanics
{
template <int DisplacementDim>
std::unique_ptr<Process> createRichardsMechanicsProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "RICHARDS_MECHANICS");
    DBUG("Create RichardsMechanicsProcess.");

    auto const staggered_scheme =
        //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(staggered_scheme && (*staggered_scheme == "staggered"));

    // Process variable.

    //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    ProcessVariable* variable_p;
    ProcessVariable* variable_u;
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    if (use_monolithic_scheme)  // monolithic scheme.
    {
        auto per_process_variables = findProcessVariables(
            variables, pv_config,
            {//! \ogs_file_param_special{prj__processes__process__RICHARDS_MECHANICS__process_variables__pressure}
             "pressure",
             //! \ogs_file_param_special{prj__processes__process__RICHARDS_MECHANICS__process_variables__displacement}
             "displacement"});
        variable_p = &per_process_variables[0].get();
        variable_u = &per_process_variables[1].get();
        process_variables.push_back(std::move(per_process_variables));
    }
    else  // staggered scheme.
    {
        using namespace std::string_literals;
        for (auto const& variable_name : {"pressure"s, "displacement"s})
        {
            auto per_process_variables =
                findProcessVariables(variables, pv_config, {variable_name});
            process_variables.push_back(std::move(per_process_variables));
        }
        variable_p = &process_variables[0][0].get();
        variable_u = &process_variables[1][0].get();
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

    DBUG("Associate pressure with process variable \'%s\'.",
         variable_p->getName().c_str());
    if (variable_p->getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Pressure process variable '%s' is not a scalar variable but has "
            "%d components.",
            variable_p->getName().c_str(),
            variable_p->getNumberOfComponents());
    }

    // Constitutive relation.
    auto solid_material =
        MaterialLib::Solids::createConstitutiveRelation<DisplacementDim>(
            parameters, config);

    // Intrinsic permeability
    auto& intrinsic_permeability = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__RICHARDS_MECHANICS__intrinsic_permeability}
        "intrinsic_permeability", parameters, 1);

    DBUG("Use \'%s\' as intrinsic conductivity parameter.",
         intrinsic_permeability.name.c_str());

    // Fluid bulk modulus
    auto& fluid_bulk_modulus = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__RICHARDS_MECHANICS__fluid_bulk_modulus}
        "fluid_bulk_modulus", parameters, 1);
    DBUG("Use \'%s\' as fluid bulk modulus parameter.",
         fluid_bulk_modulus.name.c_str());

    // Biot coefficient
    auto& biot_coefficient = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__RICHARDS_MECHANICS__biot_coefficient}
        "biot_coefficient", parameters, 1);
    DBUG("Use \'%s\' as Biot coefficient parameter.",
         biot_coefficient.name.c_str());

    // Solid density
    auto& solid_density = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__RICHARDS_MECHANICS__solid_density}
        "solid_density", parameters, 1);
    DBUG("Use \'%s\' as solid density parameter.", solid_density.name.c_str());

    // Solid bulk modulus
    auto& solid_bulk_modulus = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__RICHARDS_MECHANICS__solid_bulk_modulus}
        "solid_bulk_modulus", parameters, 1);
    DBUG("Use \'%s\' as solid bulk modulus parameter.",
         solid_bulk_modulus.name.c_str());

    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__specific_body_force}
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

    auto& temperature = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__RICHARDS_MECHANICS__temperature}
        "temperature", parameters, 1);

    auto const& flow_material_config =
        //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__material_property}
        config.getConfigSubtree("material_property");

    boost::optional<MeshLib::PropertyVector<int> const&> material_ids;
    if (mesh.getProperties().existsPropertyVector<int>("MaterialIDs"))
    {
        INFO(
            "MaterialIDs vector found; the Richards flow is in heterogeneous "
            "porous media.");
        material_ids =
            *mesh.getProperties().getPropertyVector<int>("MaterialIDs");
    }
    else
    {
        INFO(
            "MaterialIDs vector not found; the Richards flow is in homogeneous "
            "porous media.");
    }
    auto flow_material =
        ProcessLib::RichardsFlow::createRichardsFlowMaterialProperties(
            flow_material_config, material_ids, parameters);

    RichardsMechanicsProcessData<DisplacementDim> process_data{
        std::move(flow_material),
        std::move(solid_material),
        intrinsic_permeability,
        fluid_bulk_modulus,
        biot_coefficient,
        solid_density,
        solid_bulk_modulus,
        temperature,
        specific_body_force};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"RichardsMechanics_displacement"});

    ProcessLib::createSecondaryVariables(config, secondary_variables,
                                         named_function_caller);

    return std::make_unique<RichardsMechanicsProcess<DisplacementDim>>(
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller),
        use_monolithic_scheme);
}

template std::unique_ptr<Process> createRichardsMechanicsProcess<2>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<Process> createRichardsMechanicsProcess<3>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace RichardsMechanics
}  // namespace ProcessLib
