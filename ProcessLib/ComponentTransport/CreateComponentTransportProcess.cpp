/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateComponentTransportProcess.h"

#include "MaterialLib/Fluid/FluidProperties/CreateFluidProperties.h"
#include "MaterialLib/PorousMedium/CreatePorousMediaProperties.h"

#include "ProcessLib/Parameter/ConstantParameter.h"
#include "ProcessLib/Utils/ParseSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "ComponentTransportProcess.h"
#include "ComponentTransportProcessData.h"

namespace ProcessLib
{
namespace ComponentTransport
{
std::unique_ptr<Process> createComponentTransportProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "ComponentTransport");

    DBUG("Create ComponentTransportProcess.");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__ComponentTransport__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    auto process_variables = findProcessVariables(
        variables, pv_config,
        {
        //! \ogs_file_param_special{prj__processes__process__ComponentTransport__process_variables__concentration}
        "concentration",
        //! \ogs_file_param_special{prj__processes__process__ComponentTransport__process_variables__pressure}
        "pressure"});

    MaterialLib::PorousMedium::PorousMediaProperties porous_media_properties{
        MaterialLib::PorousMedium::createPorousMediaProperties(
            mesh, config)};

    //! \ogs_file_param{prj__processes__process__ComponentTransport__fluid}
    auto const& fluid_config = config.getConfigSubtree("fluid");

    auto fluid_properties =
        MaterialLib::Fluid::createFluidProperties(fluid_config);

    // Parameter for the density of the fluid.
    auto& fluid_reference_density= findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__ComponentTransport__fluid_reference_density}
        "fluid_reference_density", parameters, 1);
    DBUG("Use \'%s\' as fluid_reference_density parameter.",
         fluid_reference_density.name.c_str());

    // Parameter for the longitudinal molecular diffusion coefficient.
    auto const& molecular_diffusion_coefficient = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__ComponentTransport__molecular_diffusion_coefficient}
        "molecular_diffusion_coefficient", parameters, 1);
    DBUG("Use \'%s\' as molecular diffusion coefficient parameter.",
         molecular_diffusion_coefficient.name.c_str());

    // Parameter for the longitudinal solute dispersivity.
    auto const& solute_dispersivity_longitudinal = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__ComponentTransport__solute_dispersivity_longitudinal}
        "solute_dispersivity_longitudinal", parameters, 1);
    DBUG("Use \'%s\' as longitudinal solute dispersivity parameter.",
         solute_dispersivity_longitudinal.name.c_str());

    // Parameter for the transverse solute dispersivity.
    auto const& solute_dispersivity_transverse = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__ComponentTransport__solute_dispersivity_transverse}
        "solute_dispersivity_transverse", parameters, 1);
    DBUG("Use \'%s\' as transverse solute dispersivity parameter.",
         solute_dispersivity_transverse.name.c_str());

    // Parameter for the retardation factor.
    auto const& retardation_factor =
        findParameter<double>(config,
        //! \ogs_file_param_special{prj__processes__process__ComponentTransport__retardation_factor}
        "retardation_factor", parameters, 1);

    // Parameter for the decay rate.
    auto const& decay_rate =
        findParameter<double>(config,
        //! \ogs_file_param_special{prj__processes__process__ComponentTransport__decay_rate}
        "decay_rate", parameters, 1);

    // Specific body force parameter.
    Eigen::VectorXd specific_body_force;
    std::vector<double> const b =
        //! \ogs_file_param{prj__processes__process__ComponentTransport__specific_body_force}
        config.getConfigParameter<std::vector<double>>("specific_body_force");
    assert(b.size() > 0 && b.size() < 4);
    if (b.size() < mesh.getDimension())
        OGS_FATAL(
            "specific body force (gravity vector) has %d components, mesh "
            "dimension is %d",
            b.size(), mesh.getDimension());
    bool const has_gravity = MathLib::toVector(b).norm() > 0;
    if (has_gravity)
    {
        specific_body_force.resize(b.size());
        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    ComponentTransportProcessData process_data{
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

    NumLib::NamedFunctionCaller named_function_caller(
        {"ComponentTransport_concentration_pressure"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    return std::make_unique<ComponentTransportProcess>(
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller));
}

}  // namespace ComponentTransport
}  // namespace ProcessLib
