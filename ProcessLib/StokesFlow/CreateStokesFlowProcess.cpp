/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateStokesFlowProcess.h"
#include "StokesFlowProcess.h"

#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
namespace StokesFlow
{
namespace
{
void checkMPLProperties(
    MeshLib::Mesh const& mesh,
    MaterialPropertyLib::MaterialSpatialDistributionMap const& media_map)
{
    std::array const required_properties_liquid_phase = {
        MaterialPropertyLib::PropertyType::viscosity};

    for (auto const& element : mesh.getElements())
    {
        auto const element_id = element->getID();
        auto const& medium = *media_map.getMedium(element_id);

        // check if liquid phase definition and the corresponding properties
        // exist
        auto const& liquid_phase = medium.phase("AqueousLiquid");
        checkRequiredProperties(liquid_phase, required_properties_liquid_phase);
    }
}
}  // namespace

template <int GlobalDim>
std::unique_ptr<Process> createStokesFlowProcess(
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
    config.checkConfigParameter("type", "StokesFlow");

    DBUG("Create StokesFlowProcess.");

    auto const coupling_scheme =
        //! \ogs_file_param{prj__processes__process__StokesFlow__coupling_scheme}
        config.getConfigParameter<std::string>("coupling_scheme",
                                               "monolithic_scheme");
    const bool use_monolithic_scheme = (coupling_scheme != "staggered");

    // Process variable.
    //! \ogs_file_param{prj__processes__process__StokesFlow__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;

    // Collect all process variables in a vector before allocation
    // liquid_velocity first, pressure then
    auto const collected_process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__StokesFlow__process_variables__liquid_velocity}
         "liquid_velocity",
         //! \ogs_file_param_special{prj__processes__process__StokesFlow__process_variables__pressure}
         "pressure"});

    // Check number of components for each process variable
    auto const variable_v = collected_process_variables[0];
    if (variable_v.get().getNumberOfGlobalComponents() != GlobalDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '{:s}' is different "
            "from the global dimension: got {:d}, expected {:d}",
            variable_v.get().getName(),
            variable_v.get().getNumberOfGlobalComponents(),
            GlobalDim);
    }

    auto const variable_p = collected_process_variables[1];
    if (variable_p.get().getNumberOfGlobalComponents() != 1)
    {
        OGS_FATAL(
            "Pressure process variable '{:s}' is not a scalar variable but has "
            "{:d} components.",
            variable_p.get().getName(),
            variable_p.get().getNumberOfGlobalComponents());
    }

    // Allocate the collected process variables into a two-dimensional vector,
    // depending on what coupling scheme is adopted
    if (use_monolithic_scheme)  // monolithic scheme.
    {
        process_variables.push_back(std::move(collected_process_variables));
    }
    else  // staggered scheme.
    {
        OGS_FATAL(
            "The staggered coupling scheme for StokesFlowProcess is not "
            "implemented.");
    }

    // Specific body force
    Eigen::VectorXd specific_body_force = Eigen::VectorXd::Zero(GlobalDim);
    auto const b =
        //! \ogs_file_param{prj__processes__process__StokesFlow__specific_body_force}
        config.getConfigParameter<std::vector<double>>("specific_body_force");
    if (b.size() != GlobalDim)
    {
        OGS_FATAL(
            "The size of the specific body force vector does not match the "
            "global dimension. Vector size is {:d}, global "
            "dimension is {:d}",
            b.size(), GlobalDim);
    }
    std::copy_n(b.data(), b.size(), specific_body_force.data());

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    DBUG("Check the media properties of StokesFlow process ...");
    checkMPLProperties(mesh, *media_map);
    DBUG("Media properties verified.");

    StokesFlowProcessData process_data{std::move(media_map),
                                       std::move(specific_body_force)};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<StokesFlowProcess<GlobalDim>>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        use_monolithic_scheme);
}

template std::unique_ptr<Process> createStokesFlowProcess<2>(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media);

}  // namespace StokesFlow
}  // namespace ProcessLib
