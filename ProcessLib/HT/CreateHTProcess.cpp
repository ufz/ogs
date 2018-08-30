/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateHTProcess.h"

#include "MaterialLib/Fluid/FluidProperties/CreateFluidProperties.h"
#include "MaterialLib/PorousMedium/CreatePorousMediaProperties.h"

#include "MeshLib/IO/readMeshFromFile.h"

#include "ProcessLib/CalculateSurfaceFlux/Balance.h"
#include "ProcessLib/CalculateSurfaceFlux/ParseCalculateSurfaceFluxData.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Parameter/ConstantParameter.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "HTProcess.h"
#include "HTMaterialProperties.h"
#include "HTLocalAssemblerInterface.h"

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
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
    std::string const& output_directory)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "HT");

    DBUG("Create HTProcess.");

    auto const staggered_scheme =
        //! \ogs_file_param{prj__processes__process__HT__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(staggered_scheme && (*staggered_scheme == "staggered"));

    // Process variable.

    //! \ogs_file_param{prj__processes__process__HT__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    if (use_monolithic_scheme)  // monolithic scheme.
    {
        auto per_process_variables = findProcessVariables(
            variables, pv_config,
            {//! \ogs_file_param_special{prj__processes__process__HT__process_variables__temperature}
             "temperature",
             //! \ogs_file_param_special{prj__processes__process__HT__process_variables__pressure}
             "pressure"});
        process_variables.push_back(std::move(per_process_variables));
    }
    else  // staggered scheme.
    {
        using namespace std::string_literals;
        for (auto const& variable_name : {"temperature"s, "pressure"s})
        {
            auto per_process_variables =
                findProcessVariables(variables, pv_config, {variable_name});
            process_variables.push_back(std::move(per_process_variables));
        }
    }

    MaterialLib::PorousMedium::PorousMediaProperties porous_media_properties{
        MaterialLib::PorousMedium::createPorousMediaProperties(mesh, config,
                                                               parameters)};

    //! \ogs_file_param{prj__processes__process__HT__fluid}
    auto const& fluid_config = config.getConfigSubtree("fluid");
    auto fluid_properties =
        MaterialLib::Fluid::createFluidProperties(fluid_config);

    // Parameter for the density of the solid.
    auto& density_solid = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HT__density_solid}
        "density_solid", parameters, 1);
    DBUG("Use \'%s\' as density_solid parameter.", density_solid.name.c_str());

    // Parameter for the specific heat capacity of the solid.
    auto& specific_heat_capacity_solid = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HT__specific_heat_capacity_solid}
        "specific_heat_capacity_solid", parameters, 1);
    DBUG("Use \'%s\' as specific_heat_capacity_solid parameter.",
         specific_heat_capacity_solid.name.c_str());

    // Parameter for the thermal conductivity of the solid (only one scalar per
    // element, i.e., the isotropic case is handled at the moment)

    ConstantParameter<double> default_thermal_dispersivity_longitudinal(
        "default thermal dispersivity longitudinal", 0.);
    ConstantParameter<double> default_thermal_dispersivity_transversal(
        "default thermal dispersivity transversal", 0.);

    Parameter<double>* thermal_dispersivity_longitudinal =
        &default_thermal_dispersivity_longitudinal;
    Parameter<double>* thermal_dispersivity_transversal =
        &default_thermal_dispersivity_transversal;
    auto const dispersion_config =
        //! \ogs_file_param{prj__processes__process__HT__thermal_dispersivity}
        config.getConfigSubtreeOptional("thermal_dispersivity");
    bool const has_fluid_thermal_dispersivity =
        static_cast<bool>(dispersion_config);
    if (dispersion_config)
    {
        thermal_dispersivity_longitudinal = &findParameter<double>(
            *dispersion_config,
            //! \ogs_file_param_special{prj__processes__process__HT__thermal_dispersivity__longitudinal}
            "longitudinal", parameters, 1);
        DBUG("Use \'%s\' as thermal_dispersivity_longitudinal parameter.",
             thermal_dispersivity_longitudinal->name.c_str());

        // Parameter for the thermal conductivity of the solid (only one scalar
        // per
        // element, i.e., the isotropic case is handled at the moment)
        thermal_dispersivity_transversal = &findParameter<double>(
            *dispersion_config,
            //! \ogs_file_param_special{prj__processes__process__HT__thermal_dispersivity__transversal}
            "transversal", parameters, 1);
        DBUG("Use \'%s\' as thermal_dispersivity_transversal parameter.",
             thermal_dispersivity_transversal->name.c_str());
    }

    // Parameter for the thermal conductivity of the solid (only one scalar per
    // element, i.e., the isotropic case is handled at the moment)
    auto& thermal_conductivity_solid = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HT__thermal_conductivity_solid}
        "thermal_conductivity_solid", parameters, 1);
    DBUG("Use \'%s\' as thermal_conductivity_solid parameter.",
         thermal_conductivity_solid.name.c_str());

    // Parameter for the thermal conductivity of the fluid.
    auto& thermal_conductivity_fluid = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__HT__thermal_conductivity_fluid}
        "thermal_conductivity_fluid", parameters, 1);
    DBUG("Use \'%s\' as thermal_conductivity_fluid parameter.",
         thermal_conductivity_fluid.name.c_str());

    // Specific body force parameter.
    Eigen::VectorXd specific_body_force;
    std::vector<double> const b =
        //! \ogs_file_param{prj__processes__process__HT__specific_body_force}
        config.getConfigParameter<std::vector<double>>("specific_body_force");
    assert(!b.empty() && b.size() < 4);
    if (b.size() < mesh.getDimension())
    {
        OGS_FATAL(
            "specific body force (gravity vector) has %d components, mesh "
            "dimension is %d",
            b.size(), mesh.getDimension());
    }
    bool const has_gravity = MathLib::toVector(b).norm() > 0;
    if (has_gravity)
    {
        specific_body_force.resize(b.size());
        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    ConstantParameter<double> default_solid_thermal_expansion(
        "default solid thermal expansion", 0.);
    ConstantParameter<double> default_biot_constant("default_biot constant",
                                                    0.);
    Parameter<double>* solid_thermal_expansion =
        &default_solid_thermal_expansion;
    Parameter<double>* biot_constant = &default_biot_constant;

    auto const solid_config =
        //! \ogs_file_param{prj__processes__process__HT__solid_thermal_expansion}
        config.getConfigSubtreeOptional("solid_thermal_expansion");
    const bool has_fluid_thermal_expansion = static_cast<bool>(solid_config);
    if (solid_config)
    {
        solid_thermal_expansion = &findParameter<double>(
            //! \ogs_file_param_special{prj__processes__process__HT__solid_thermal_expansion__thermal_expansion}
            *solid_config, "thermal_expansion", parameters, 1);
        DBUG("Use \'%s\' as solid thermal expansion.",
             solid_thermal_expansion->name.c_str());
        biot_constant = &findParameter<double>(
            //! \ogs_file_param_special{prj__processes__process__HT__solid_thermal_expansion__biot_constant}
            *solid_config, "biot_constant", parameters, 1);
        DBUG("Use \'%s\' as Biot's constant.", biot_constant->name.c_str());
    }

    // for the balance
    std::unique_ptr<ProcessLib::Balance> balance;
    std::string mesh_name;  // surface mesh the balance will computed on
    std::string balance_pv_name;
    std::string balance_out_fname;
    ProcessLib::parseCalculateSurfaceFluxData(
        config, mesh_name, balance_pv_name, balance_out_fname);

    if (!mesh_name.empty())  // balance is optional
    {
        balance_out_fname =
            BaseLib::copyPathToFileName(balance_out_fname, output_directory);

        balance.reset(new ProcessLib::Balance(std::move(mesh_name), meshes,
                                              std::move(balance_pv_name),
                                              std::move(balance_out_fname)));

        // Surface mesh and bulk mesh must have equal axial symmetry flags!
        balance->surface_mesh.setAxiallySymmetric(mesh.isAxiallySymmetric());
    }

    std::unique_ptr<HTMaterialProperties> material_properties =
        std::make_unique<HTMaterialProperties>(
            std::move(porous_media_properties),
            density_solid,
            std::move(fluid_properties),
            has_fluid_thermal_dispersivity,
            *thermal_dispersivity_longitudinal,
            *thermal_dispersivity_transversal,
            specific_heat_capacity_solid,
            thermal_conductivity_solid,
            thermal_conductivity_fluid,
            has_fluid_thermal_expansion,
            *solid_thermal_expansion,
            *biot_constant,
            specific_body_force,
            has_gravity);

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"HT_temperature_pressure"});

    ProcessLib::createSecondaryVariables(config, secondary_variables,
                                         named_function_caller);

    return std::make_unique<HTProcess>(
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(material_properties),
        std::move(secondary_variables), std::move(named_function_caller),
        use_monolithic_scheme, std::move(balance));
}

}  // namespace HT
}  // namespace ProcessLib
