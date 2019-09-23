/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateComponentTransportProcess.h"

#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/SurfaceFlux/SurfaceFluxData.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "ComponentTransportProcess.h"
#include "ComponentTransportProcessData.h"
namespace ProcessLib
{
namespace ComponentTransport
{
std::unique_ptr<Process> createComponentTransportProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
    std::string const& output_directory,
    std::map<int, std::unique_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "ComponentTransport");

    DBUG("Create ComponentTransportProcess.");
    auto const staggered_scheme =
        //! \ogs_file_param{prj__processes__process__ComponentTransport__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(staggered_scheme && (*staggered_scheme == "staggered"));

    // Process variable.

    //! \ogs_file_param{prj__processes__process__ComponentTransport__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;

    // Collect all process variables in a vector before allocation
    // pressure first, concentration then
    auto const collected_process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__ComponentTransport__process_variables__pressure}
         "pressure",
         //! \ogs_file_param_special{prj__processes__process__ComponentTransport__process_variables__concentration}
         "concentration"});

    // Check number of components for each process variable
    auto it = std::find_if(
        collected_process_variables.cbegin(),
        collected_process_variables.cend(),
        [](std::reference_wrapper<ProcessLib::ProcessVariable> const& pv) {
            return pv.get().getNumberOfComponents() != 1;
        });

    if (it != collected_process_variables.end())
    {
        OGS_FATAL(
            "Number of components for process variable '%s' should be 1 rather "
            "than %d.",
            it->get().getName().c_str(),
            it->get().getNumberOfComponents());
    }

    std::vector<std::pair<int, std::string>> process_id_to_component_name_map;
    // Allocate the collected process variables into a two-dimensional vector,
    // depending on what scheme is adopted
    if (use_monolithic_scheme)  // monolithic scheme.
    {
        process_variables.push_back(std::move(collected_process_variables));
    }
    else  // staggered scheme.
    {
        std::vector<std::reference_wrapper<ProcessLib::ProcessVariable>>
            per_process_variable;

        for (auto& pv : collected_process_variables)
        {
            per_process_variable.emplace_back(pv);
            process_variables.push_back(std::move(per_process_variable));
        }

        auto variable_id = 0;
        for (unsigned process_id = 1; process_id < process_variables.size();
             process_id++)
        {
            auto const& transport_process_variable =
                process_variables[process_id][variable_id].get().getName();
            process_id_to_component_name_map.emplace_back(
                process_id, transport_process_variable);
        }
    }

    // Specific body force parameter.
    Eigen::VectorXd specific_body_force;
    std::vector<double> const b =
        //! \ogs_file_param{prj__processes__process__ComponentTransport__specific_body_force}
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

    bool const non_advective_form =
        //! \ogs_file_param{prj__processes__process__ComponentTransport__non_advective_form}
        config.getConfigParameter<bool>("non_advective_form", false);

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    ComponentTransportProcessData process_data{
        std::move(media_map),
        specific_body_force,
        has_gravity,
        non_advective_form};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"ComponentTransport_concentration_pressure"});

    ProcessLib::createSecondaryVariables(config, secondary_variables,
                                         named_function_caller);

    std::unique_ptr<ProcessLib::SurfaceFluxData> surfaceflux;
    auto surfaceflux_config =
        //! \ogs_file_param{prj__processes__process__calculatesurfaceflux}
        config.getConfigSubtreeOptional("calculatesurfaceflux");
    if (surfaceflux_config)
    {
        surfaceflux = ProcessLib::SurfaceFluxData::createSurfaceFluxData(
            *surfaceflux_config, meshes, output_directory);
    }

    return std::make_unique<ComponentTransportProcess>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        std::move(named_function_caller), use_monolithic_scheme,
        std::move(surfaceflux), std::move(process_id_to_component_name_map));
}

}  // namespace ComponentTransport
}  // namespace ProcessLib
