/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateHTProcess.h"

#include "MaterialLib/MPL/CheckMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "ParameterLib/ConstantParameter.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/SurfaceFlux/SurfaceFluxData.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "HTProcess.h"
#include "HTProcessData.h"
#include "HTLocalAssemblerInterface.h"

namespace ProcessLib
{
namespace HT
{
void checkMPLProperties(
    MeshLib::Mesh const& mesh,
    MaterialPropertyLib::MaterialSpatialDistributionMap const& media_map)
{
    std::array const required_property_medium = {
        MaterialPropertyLib::PropertyType::porosity,
        MaterialPropertyLib::PropertyType::permeability};

    std::array const required_property_liquid_phase = {
        MaterialPropertyLib::PropertyType::viscosity,
        MaterialPropertyLib::PropertyType::density,
        MaterialPropertyLib::PropertyType::specific_heat_capacity};

    std::array const required_property_solid_phase = {
        MaterialPropertyLib::PropertyType::specific_heat_capacity,
        MaterialPropertyLib::PropertyType::density,
        MaterialPropertyLib::PropertyType::storage};

    MaterialPropertyLib::checkMaterialSpatialDistributionMap(
        mesh, media_map, required_property_medium, required_property_solid_phase,
        required_property_liquid_phase);
}

std::unique_ptr<Process> createHTProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
    std::string const& output_directory,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
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

    // Process IDs, which are set according to the appearance order of the
    // process variables.
    int const heat_transport_process_id = 0;
    int hydraulic_process_id = 0;

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
        hydraulic_process_id = 1;
    }

    // Specific body force parameter.
    Eigen::VectorXd specific_body_force;
    std::vector<double> const b =
        //! \ogs_file_param{prj__processes__process__HT__specific_body_force}
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

    ParameterLib::ConstantParameter<double> default_solid_thermal_expansion(
        "default solid thermal expansion", 0.);
    ParameterLib::ConstantParameter<double> default_biot_constant(
        "default_biot constant", 0.);
    ParameterLib::Parameter<double>* solid_thermal_expansion =
        &default_solid_thermal_expansion;
    ParameterLib::Parameter<double>* biot_constant = &default_biot_constant;

    auto const solid_config =
        //! \ogs_file_param{prj__processes__process__HT__solid_thermal_expansion}
        config.getConfigSubtreeOptional("solid_thermal_expansion");
    const bool has_fluid_thermal_expansion = static_cast<bool>(solid_config);
    if (solid_config)
    {
        solid_thermal_expansion = &ParameterLib::findParameter<double>(
            //! \ogs_file_param_special{prj__processes__process__HT__solid_thermal_expansion__thermal_expansion}
            *solid_config, "thermal_expansion", parameters, 1, &mesh);
        DBUG("Use '{:s}' as solid thermal expansion.",
             solid_thermal_expansion->name);
        biot_constant = &ParameterLib::findParameter<double>(
            //! \ogs_file_param_special{prj__processes__process__HT__solid_thermal_expansion__biot_constant}
            *solid_config, "biot_constant", parameters, 1, &mesh);
        DBUG("Use '{:s}' as Biot's constant.", biot_constant->name);
    }

    std::unique_ptr<ProcessLib::SurfaceFluxData> surfaceflux;
    auto calculatesurfaceflux_config =
        //! \ogs_file_param{prj__processes__process__calculatesurfaceflux}
        config.getConfigSubtreeOptional("calculatesurfaceflux");
    if (calculatesurfaceflux_config)
    {
        surfaceflux = ProcessLib::SurfaceFluxData::createSurfaceFluxData(
            *calculatesurfaceflux_config, meshes);
    }

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    DBUG("Check the media properties of HT process ...");
    checkMPLProperties(mesh, *media_map);
    DBUG("Media properties verified.");

    HTProcessData process_data{
        std::move(media_map),      has_fluid_thermal_expansion,
        *solid_thermal_expansion,  *biot_constant,
        specific_body_force,       has_gravity,
        heat_transport_process_id, hydraulic_process_id};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<HTProcess>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        use_monolithic_scheme, std::move(surfaceflux));
}

}  // namespace HT
}  // namespace ProcessLib
