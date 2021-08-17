/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateComponentTransportProcess.h"

#include "ChemistryLib/ChemicalSolverInterface.h"
#include "ComponentTransportProcess.h"
#include "ComponentTransportProcessData.h"
#include "CreateLookupTable.h"
#include "LookupTable.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/SurfaceFlux/SurfaceFluxData.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
namespace ComponentTransport
{
void checkMPLProperties(
    MeshLib::Mesh const& mesh,
    MaterialPropertyLib::MaterialSpatialDistributionMap const& media_map)
{
    std::array const required_properties_medium = {
        MaterialPropertyLib::PropertyType::porosity,
        MaterialPropertyLib::PropertyType::transversal_dispersivity,
        MaterialPropertyLib::PropertyType::longitudinal_dispersivity,
        MaterialPropertyLib::PropertyType::permeability};

    std::array const required_properties_liquid_phase = {
        MaterialPropertyLib::PropertyType::density,
        MaterialPropertyLib::PropertyType::viscosity};

    std::array const required_properties_components = {
        MaterialPropertyLib::PropertyType::retardation_factor,
        MaterialPropertyLib::PropertyType::decay_rate,
        MaterialPropertyLib::PropertyType::pore_diffusion};

    for (auto const& element : mesh.getElements())
    {
        auto const element_id = element->getID();

        auto const& medium = *media_map.getMedium(element_id);
        checkRequiredProperties(medium, required_properties_medium);

        // check if liquid phase definition and the corresponding properties
        // exist
        auto const& liquid_phase = medium.phase("AqueousLiquid");
        checkRequiredProperties(liquid_phase, required_properties_liquid_phase);

        // check if components and the corresponding properties exist
        auto const number_of_components = liquid_phase.numberOfComponents();
        for (std::size_t component_id = 0; component_id < number_of_components;
             ++component_id)
        {
            if (!liquid_phase.hasComponent(component_id))
            {
                OGS_FATAL(
                    "The component {:d} in the AqueousLiquid phase isn't "
                    "specified.",
                    component_id);
            }
            auto const& component = liquid_phase.component(component_id);
            checkRequiredProperties(component, required_properties_components);
        }
    }
}

std::unique_ptr<Process> createComponentTransportProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media,
    std::unique_ptr<ChemistryLib::ChemicalSolverInterface>&&
        chemical_solver_interface)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "ComponentTransport");

    DBUG("Create ComponentTransportProcess.");

    auto const coupling_scheme =
        //! \ogs_file_param{prj__processes__process__ComponentTransport__coupling_scheme}
        config.getConfigParameter<std::string>("coupling_scheme",
                                               "monolithic_scheme");
    const bool use_monolithic_scheme = (coupling_scheme != "staggered");

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
        [](std::reference_wrapper<ProcessLib::ProcessVariable> const& pv)
        { return pv.get().getNumberOfGlobalComponents() != 1; });

    if (it != collected_process_variables.end())
    {
        OGS_FATAL(
            "Number of components for process variable '{:s}' should be 1 "
            "rather "
            "than {:d}.",
            it->get().getName(),
            it->get().getNumberOfGlobalComponents());
    }
    int const hydraulic_process_id = 0;
    int const first_transport_process_id = use_monolithic_scheme ? 0 : 1;

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

        if (!chemical_solver_interface)
        {
            for (auto& pv : collected_process_variables)
            {
                per_process_variable.emplace_back(pv);
                process_variables.push_back(std::move(per_process_variable));
            }
        }
        else
        {
            auto sort_by_component =
                [&per_process_variable,
                 collected_process_variables](auto const& c_name)
            {
                auto pv = std::find_if(collected_process_variables.begin(),
                                       collected_process_variables.end(),
                                       [&c_name](auto const& v) -> bool
                                       { return v.get().getName() == c_name; });

                if (pv == collected_process_variables.end())
                {
                    OGS_FATAL(
                        "Component {:s} given in "
                        "<chemical_system>/<solution>/"
                        "<components> is not found in specified "
                        "coupled processes (see "
                        "<process>/<process_variables>/"
                        "<concentration>).",
                        c_name);
                }

                per_process_variable.emplace_back(*pv);
                return std::move(per_process_variable);
            };

            auto const components =
                chemical_solver_interface->getComponentList();
            // pressure
            per_process_variable.emplace_back(collected_process_variables[0]);
            process_variables.push_back(std::move(per_process_variable));
            // concentration
            assert(components.size() + 1 == collected_process_variables.size());
            std::transform(components.begin(), components.end(),
                           std::back_inserter(process_variables),
                           sort_by_component);
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

    bool const non_advective_form =
        //! \ogs_file_param{prj__processes__process__ComponentTransport__non_advective_form}
        config.getConfigParameter<bool>("non_advective_form", false);

    bool chemically_induced_porosity_change =
        //! \ogs_file_param{prj__processes__process__ComponentTransport__chemically_induced_porosity_change}
        config.getConfigParameter<bool>("chemically_induced_porosity_change",
                                        false);

    auto const temperature = ParameterLib::findOptionalTagParameter<double>(
        //! \ogs_file_param_special{prj__processes__process__ComponentTransport__temperature_field}
        config, "temperature_field", parameters, 1, &mesh);

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    auto lookup_table = ComponentTransport::createLookupTable(
        //! \ogs_file_param_special{prj__processes__process__ComponentTransport__tabular_file}
        config.getConfigParameterOptional<std::string>("tabular_file"),
        process_variables);

    DBUG("Check the media properties of ComponentTransport process ...");
    checkMPLProperties(mesh, *media_map);
    DBUG("Media properties verified.");

    ComponentTransportProcessData process_data{
        std::move(media_map),
        specific_body_force,
        has_gravity,
        non_advective_form,
        temperature,
        chemically_induced_porosity_change,
        chemical_solver_interface.get(),
        std::move(lookup_table),
        hydraulic_process_id,
        first_transport_process_id};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    std::unique_ptr<ProcessLib::SurfaceFluxData> surfaceflux;
    auto surfaceflux_config =
        //! \ogs_file_param{prj__processes__process__calculatesurfaceflux}
        config.getConfigSubtreeOptional("calculatesurfaceflux");
    if (surfaceflux_config)
    {
        surfaceflux = ProcessLib::SurfaceFluxData::createSurfaceFluxData(
            *surfaceflux_config, meshes);
    }

    return std::make_unique<ComponentTransportProcess>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        use_monolithic_scheme, std::move(surfaceflux),
        std::move(chemical_solver_interface));
}

}  // namespace ComponentTransport
}  // namespace ProcessLib
