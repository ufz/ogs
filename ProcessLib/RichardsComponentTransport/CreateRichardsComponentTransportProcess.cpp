/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateRichardsComponentTransportProcess.h"

#include "MaterialLib/Component/CreateComponentProperties.h"
#include "MaterialLib/Fluid/FluidProperties/CreateFluidProperties.h"

#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Parameter/ConstantParameter.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "CreatePorousMediaProperties.h"
#include "RichardsComponentTransportProcess.h"
#include "RichardsComponentTransportProcessData.h"

namespace ProcessLib
{
namespace RichardsComponentTransport
{
std::unique_ptr<Process> createRichardsComponentTransportProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
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
            {
            //! \ogs_file_param_special{prj__processes__process__RichardsComponentTransport__process_variables__concentration}
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
    auto& fluid_reference_density= findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__RichardsComponentTransport__fluid_reference_density}
        "fluid_reference_density", parameters, 1);
    DBUG("Use \'%s\' as fluid_reference_density parameter.",
         fluid_reference_density.name.c_str());

    auto component_properties =
        MaterialLib::Component::createComponentProperties(config, parameters);

    {
        auto variable = std::find_if(variables.cbegin(), variables.cend(),
                                     [](ProcessLib::ProcessVariable const& v) {
                                         return v.getName() == "concentration";
                                     });

        if ((unsigned) variable->getNumberOfComponents() != component_properties.size())
            OGS_FATAL(
                "'Concentration' has %d component(s), while %d types of "
                "component(s) are prescribed.",
                variable->getNumberOfComponents(),
                component_properties.size());
    }

    // Parameter for the retardation factor.
    auto const& retardation_factor = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__RichardsComponentTransport__retardation_factor}
        "retardation_factor", parameters, 1);

    // Parameter for the decay rate.
    auto const& decay_rate = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__RichardsComponentTransport__decay_rate}
        "decay_rate", parameters, 1);

    // Specific body force parameter.
    Eigen::VectorXd specific_body_force;
    std::vector<double> const b =
        //! \ogs_file_param{prj__processes__process__RichardsComponentTransport__specific_body_force}
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

    RichardsComponentTransportProcessData process_data{
        std::move(porous_media_properties),
        fluid_reference_density,
        std::move(fluid_properties),
        std::move(component_properties),
        retardation_factor,
        decay_rate,
        specific_body_force,
        has_gravity};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"RichardsComponentTransport_concentration_pressure"});

    ProcessLib::createSecondaryVariables(config, secondary_variables,
                                         named_function_caller);

    return std::make_unique<RichardsComponentTransportProcess>(
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller),
        use_monolithic_scheme);
}

}  // namespace RichardsComponentTransport
}  // namespace ProcessLib
