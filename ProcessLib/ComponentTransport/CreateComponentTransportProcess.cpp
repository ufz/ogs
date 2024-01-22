/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "MeshLib/Utils/GetElementRotationMatrices.h"
#include "MeshLib/Utils/GetSpaceDimension.h"
#include "NumLib/NumericalStability/CreateNumericalStabilization.h"
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

    for (auto const& element_id : mesh.getElements() | MeshLib::views::ids)
    {
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
    std::string const& name,
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

    /// \section processvariablescomponenttransport Process Variables

    //! \ogs_file_param{prj__processes__process__ComponentTransport__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;

    /// Primary process variables as they appear in the global component vector.
    std::vector collected_process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__ComponentTransport__process_variables__pressure}
         "pressure",
         //! \ogs_file_param_special{prj__processes__process__ComponentTransport__process_variables__concentration}
         "concentration"});

    /// Temperature is optional.
    std::vector const temperature_variable = findProcessVariables(
        variables, pv_config,
        //! \ogs_file_param_special{prj__processes__process__ComponentTransport__process_variables__temperature}
        "temperature", true /*temperature is optional*/);
    bool const isothermal = temperature_variable.empty();
    if (!isothermal)
    {
        assert(temperature_variable.size() == 1);
        collected_process_variables.insert(
            ++collected_process_variables.begin(), temperature_variable[0]);
    }

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

    // Allocate the collected process variables into a two-dimensional vector,
    // depending on what scheme is adopted
    if (use_monolithic_scheme)  // monolithic scheme.
    {
        if (!isothermal)
        {
            OGS_FATAL(
                "Currently, non-isothermal component transport process can "
                "only be simulated in staggered scheme.");
        }

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
            // temperature
            if (!isothermal)
            {
                per_process_variable.emplace_back(
                    collected_process_variables[1]);
                process_variables.push_back(std::move(per_process_variable));
            }
            // concentration
            assert(components.size() + (isothermal ? 1 : 2) ==
                   collected_process_variables.size());
            std::transform(components.begin(), components.end(),
                           std::back_inserter(process_variables),
                           sort_by_component);
        }
    }

    /// \section parameterscomponenttransport Process Parameters
    std::vector<double> const b =
        //! \ogs_file_param{prj__processes__process__ComponentTransport__specific_body_force}
        config.getConfigParameter<std::vector<double>>("specific_body_force");
    assert(!b.empty() && b.size() < 4);
    // Specific body force parameter.
    Eigen::VectorXd specific_body_force(b.size());
    int const mesh_space_dimension =
        MeshLib::getSpaceDimension(mesh.getNodes());
    if (static_cast<int>(b.size()) != mesh_space_dimension)
    {
        OGS_FATAL(
            "specific body force (gravity vector) has {:d} components, mesh "
            "dimension is {:d}",
            b.size(), mesh_space_dimension);
    }
    bool const has_gravity = MathLib::toVector(b).norm() > 0;
    if (has_gravity)
    {
        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    bool const non_advective_form =
        //! \ogs_file_param{prj__processes__process__ComponentTransport__non_advective_form}
        config.getConfigParameter<bool>("non_advective_form", false);

    bool chemically_induced_porosity_change =
        //! \ogs_file_param{prj__processes__process__ComponentTransport__chemically_induced_porosity_change}
        config.getConfigParameter<bool>("chemically_induced_porosity_change",
                                        false);

    auto const temperature_field = ParameterLib::findOptionalTagParameter<
        double>(
        //! \ogs_file_param_special{prj__processes__process__ComponentTransport__temperature_field}
        config, "temperature_field", parameters, 1, &mesh);
    if (!isothermal && temperature_field != nullptr)
    {
        OGS_FATAL("Temperature field is set for non-isothermal setup.")
    }

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    auto lookup_table = ComponentTransport::createLookupTable(
        //! \ogs_file_param{prj__processes__process__ComponentTransport__tabular_file}
        config.getConfigParameterOptional<std::string>("tabular_file"),
        process_variables);

    DBUG("Check the media properties of ComponentTransport process ...");
    checkMPLProperties(mesh, media_map);
    DBUG("Media properties verified.");

    auto stabilizer = NumLib::createNumericalStabilization(mesh, config);

    auto const* aperture_size_parameter = &ParameterLib::findParameter<double>(
        ProcessLib::Process::constant_one_parameter_name, parameters, 1);
    auto const aperture_config =
        //! \ogs_file_param{prj__processes__process__ComponentTransport__aperture_size}
        config.getConfigSubtreeOptional("aperture_size");
    if (aperture_config)
    {
        aperture_size_parameter = &ParameterLib::findParameter<double>(
            //! \ogs_file_param_special{prj__processes__process__ComponentTransport__aperture_size__parameter}
            *aperture_config, "parameter", parameters, 1);
    }

    auto const is_linear =
        //! \ogs_file_param{prj__processes__process__ComponentTransport__linear}
        config.getConfigParameter("linear", false);

    auto const ls_compute_only_upon_timestep_change =
        //! \ogs_file_param{prj__processes__process__ComponentTransport__linear_solver_compute_only_upon_timestep_change}
        config.getConfigParameter(
            "linear_solver_compute_only_upon_timestep_change", false);

    auto const rotation_matrices = MeshLib::getElementRotationMatrices(
        mesh_space_dimension, mesh.getDimension(), mesh.getElements());
    std::vector<Eigen::VectorXd> projected_specific_body_force_vectors;
    projected_specific_body_force_vectors.reserve(rotation_matrices.size());

    std::transform(rotation_matrices.begin(), rotation_matrices.end(),
                   std::back_inserter(projected_specific_body_force_vectors),
                   [&specific_body_force](const auto& R)
                   { return R * R.transpose() * specific_body_force; });

    ComponentTransportProcessData process_data{
        std::move(media_map),
        has_gravity,
        non_advective_form,
        temperature_field,
        chemically_induced_porosity_change,
        chemical_solver_interface.get(),
        std::move(lookup_table),
        std::move(stabilizer),
        projected_specific_body_force_vectors,
        mesh_space_dimension,
        *aperture_size_parameter,
        isothermal,
    };

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
        std::move(chemical_solver_interface), is_linear,
        ls_compute_only_upon_timestep_change);
}

}  // namespace ComponentTransport
}  // namespace ProcessLib
