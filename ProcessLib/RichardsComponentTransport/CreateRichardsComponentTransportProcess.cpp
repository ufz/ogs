/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateRichardsComponentTransportProcess.h"
#include "CreatePorousMediaProperties.h"

#include "MaterialLib/Fluid/FluidProperties/CreateFluidProperties.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "ParameterLib/ConstantParameter.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "RichardsComponentTransportProcess.h"
#include "RichardsComponentTransportProcessData.h"

namespace ProcessLib
{
namespace RichardsComponentTransport
{
std::unique_ptr<Process> createRichardsComponentTransportProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "RichardsComponentTransport");

    DBUG("Create RichardsComponentTransportProcess.");

    auto const staggered_scheme =
        //! \ogs_file_param{prj__processes__process__RichardsComponentTransport__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(staggered_scheme && (*staggered_scheme == "staggered"));

    // Process variable.

    //! \ogs_file_param{prj__processes__process__RichardsComponentTransport__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    if (use_monolithic_scheme)  // monolithic scheme.
    {
        auto per_process_variables = findProcessVariables(
            variables, pv_config,
            {//! \ogs_file_param_special{prj__processes__process__RichardsComponentTransport__process_variables__concentration}
             "concentration",
             //! \ogs_file_param_special{prj__processes__process__RichardsComponentTransport__process_variables__pressure}
             "pressure"});
        process_variables.push_back(std::move(per_process_variables));
    }
    else  // staggered scheme.
    {
        using namespace std::string_literals;
        for (auto const& variable_name : {"concentration"s, "pressure"s})
        {
            auto per_process_variables =
                findProcessVariables(variables, pv_config, {variable_name});
            process_variables.push_back(std::move(per_process_variables));
        }
    }

    auto const& porous_medium_configs =
        //! \ogs_file_param{prj__processes__process__RichardsComponentTransport__porous_medium}
        config.getConfigSubtree("porous_medium");
    PorousMediaProperties porous_media_properties{
        createPorousMediaProperties(mesh, porous_medium_configs, parameters)};

    //! \ogs_file_param{prj__processes__process__RichardsComponentTransport__fluid}
    auto const& fluid_config = config.getConfigSubtree("fluid");

    auto fluid_properties =
        MaterialLib::Fluid::createFluidProperties(fluid_config);

    // Parameter for the density of the fluid.
    auto const& fluid_reference_density = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__RichardsComponentTransport__fluid_reference_density}
        "fluid_reference_density", parameters, 1, &mesh);
    DBUG("Use '{:s}' as fluid_reference_density parameter.",
         fluid_reference_density.name);

    // Parameter for the longitudinal solute dispersivity.
    auto const& molecular_diffusion_coefficient = ParameterLib::findParameter<
        double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__RichardsComponentTransport__molecular_diffusion_coefficient}
        "molecular_diffusion_coefficient", parameters, 1, &mesh);
    DBUG("Use '{:s}' as molecular diffusion coefficient parameter.",
         molecular_diffusion_coefficient.name);

    // Parameter for the longitudinal solute dispersivity.
    auto const& solute_dispersivity_longitudinal = ParameterLib::findParameter<
        double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__RichardsComponentTransport__solute_dispersivity_longitudinal}
        "solute_dispersivity_longitudinal", parameters, 1, &mesh);
    DBUG("Use '{:s}' as longitudinal solute dispersivity parameter.",
         solute_dispersivity_longitudinal.name);

    // Parameter for the transverse solute dispersivity.
    auto const& solute_dispersivity_transverse = ParameterLib::findParameter<
        double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__RichardsComponentTransport__solute_dispersivity_transverse}
        "solute_dispersivity_transverse", parameters, 1, &mesh);
    DBUG("Use '{:s}' as transverse solute dispersivity parameter.",
         solute_dispersivity_transverse.name);

    // Parameter for the retardation factor.
    auto const& retardation_factor = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__RichardsComponentTransport__retardation_factor}
        "retardation_factor", parameters, 1, &mesh);

    // Parameter for the decay rate.
    auto const& decay_rate = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__RichardsComponentTransport__decay_rate}
        "decay_rate", parameters, 1, &mesh);

    // Specific body force parameter.
    Eigen::VectorXd specific_body_force;
    std::vector<double> const b =
        //! \ogs_file_param{prj__processes__process__RichardsComponentTransport__specific_body_force}
        config.getConfigParameter<std::vector<double>>("specific_body_force");
    assert(!b.empty() && b.size() < 4);
    if (b.size() < mesh.getDimension())
    {
        OGS_FATAL(
            "specific body force (gravity vector) has {:d} components, mesh "
            "dimension is {:d}",
            b.size(), mesh.getDimension());
    }
    bool const has_gravity = MathLib::toVector(b).norm() > 0;
    if (has_gravity)
    {
        specific_body_force.resize(b.size());
        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    RichardsComponentTransportProcessData process_data{
        std::move(media_map),
        std::move(porous_media_properties),
        fluid_reference_density,
        std::move(fluid_properties),
        molecular_diffusion_coefficient,
        solute_dispersivity_longitudinal,
        solute_dispersivity_transverse,
        retardation_factor,
        decay_rate,
        specific_body_force,
        has_gravity};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<RichardsComponentTransportProcess>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        use_monolithic_scheme);
}

}  // namespace RichardsComponentTransport
}  // namespace ProcessLib
