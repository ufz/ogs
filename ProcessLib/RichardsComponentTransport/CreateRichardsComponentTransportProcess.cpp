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

#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "RichardsComponentTransportProcess.h"
#include "RichardsComponentTransportProcessData.h"

namespace ProcessLib
{
namespace RichardsComponentTransport
{
namespace
{
void checkMPLProperties(
    MeshLib::Mesh const& mesh,
    MaterialPropertyLib::MaterialSpatialDistributionMap const& media_map)
{
    std::array const required_properties_medium = {
        MaterialPropertyLib::PropertyType::porosity,
        MaterialPropertyLib::PropertyType::transversal_dispersivity,
        MaterialPropertyLib::PropertyType::longitudinal_dispersivity,
        MaterialPropertyLib::PropertyType::permeability,
        MaterialPropertyLib::PropertyType::saturation,
        MaterialPropertyLib::PropertyType::storage,
        MaterialPropertyLib::PropertyType::relative_permeability};

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
}  // namespace

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

    auto const coupling_scheme =
        //! \ogs_file_param{prj__processes__process__RichardsComponentTransport__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(coupling_scheme && (*coupling_scheme == "staggered"));

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
        if (per_process_variables.size() > 2)
        {
            OGS_FATAL(
                "By now RichardsComponentTransport process only supports "
                "single component transport simulation.");
        }
        process_variables.push_back(std::move(per_process_variables));
    }
    else  // staggered scheme.
    {
        OGS_FATAL("The staggered coupling scheme is not implemented.");
    }

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

    DBUG(
        "Check the media properties of RichardsComponentTransport process ...");
    checkMPLProperties(mesh, *media_map);
    DBUG("Media properties verified.");

    RichardsComponentTransportProcessData process_data{
        std::move(media_map),
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
