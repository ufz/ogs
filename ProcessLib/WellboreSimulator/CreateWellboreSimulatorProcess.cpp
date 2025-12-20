// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateWellboreSimulatorProcess.h"

#include "MaterialLib/MPL/CheckMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "ParameterLib/ConstantParameter.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "ReservoirProperties.h"
#include "WellboreGeometry.h"
#include "WellboreSimulatorLocalAssemblerInterface.h"
#include "WellboreSimulatorProcess.h"
#include "WellboreSimulatorProcessData.h"

namespace ProcessLib
{
namespace WellboreSimulator
{
void checkMPLProperties(
    MeshLib::Mesh const& mesh,
    MaterialPropertyLib::MaterialSpatialDistributionMap const& media_map)
{
    std::array<MaterialPropertyLib::PropertyType, 0> const
        required_property_medium{};

    std::array<MaterialPropertyLib::PropertyType, 0> const
        required_property_solid_phase{};

    std::array const required_property_liquid_phase = {
        MaterialPropertyLib::PropertyType::viscosity,
        MaterialPropertyLib::PropertyType::density,
        MaterialPropertyLib::PropertyType::saturation_density,
        MaterialPropertyLib::PropertyType::saturation_enthalpy,
        MaterialPropertyLib::PropertyType::temperature,
        MaterialPropertyLib::PropertyType::enthalpy};

    std::array const required_property_gas_phase = {
        MaterialPropertyLib::PropertyType::saturation_density,
        MaterialPropertyLib::PropertyType::saturation_enthalpy,
        MaterialPropertyLib::PropertyType::saturation_temperature,
        MaterialPropertyLib::PropertyType::viscosity};

    MaterialPropertyLib::checkMaterialSpatialDistributionMap(
        mesh, media_map, required_property_medium,
        required_property_solid_phase, required_property_liquid_phase,
        required_property_gas_phase);
}

std::unique_ptr<Process> createWellboreSimulatorProcess(
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
    config.checkConfigParameter("type", "WELLBORE_SIMULATOR");

    DBUG("Create WellboreSimulatorProcess.");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__WELLBORE_SIMULATOR__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;

    auto collected_process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__WELLBORE_SIMULATOR__process_variables__pressure}
         "pressure",
         //! \ogs_file_param_special{prj__processes__process__WELLBORE_SIMULATOR__process_variables__velocity}
         "velocity",
         //! \ogs_file_param_special{prj__processes__process__WELLBORE_SIMULATOR__process_variables__enthalpy}
         "specific_enthalpy"});

    process_variables.push_back(std::move(collected_process_variables));

    // Specific body force parameter.
    Eigen::VectorXd specific_body_force;
    std::vector<double> const b =
        //! \ogs_file_param{prj__processes__process__WELLBORE_SIMULATOR__specific_body_force}
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

    WellboreGeometry const wellbore_geometry =
        //! \ogs_file_param{prj__processes__process__WELLBORE_SIMULATOR__wellbore}
        createWellboreGeometry(config.getConfigSubtree("wellbore"), parameters);

    auto const& well_ref_pressure = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__WELLBORE_SIMULATOR__wellbore_ref_pressure}
        "wellbore_ref_pressure", parameters, 1, &mesh);
    DBUG("Use '{:s}' as wellbore_ref_pressure parameter.",
         well_ref_pressure.name);

    auto const& well_ref_enthalpy = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__WELLBORE_SIMULATOR__wellbore_ref_enthalpy}
        "wellbore_ref_enthalpy", parameters, 1, &mesh);
    DBUG("Use '{:s}' as wellbore_ref_enthalpy parameter.",
         well_ref_enthalpy.name);

    auto const heat_exchange_with_formation =
        //! \ogs_file_param{prj__processes__process__WELLBORE_SIMULATOR__heat_exchange_with_formation}
        config.getConfigParameter<bool>("heat_exchange_with_formation", false);

    ReservoirProperties const reservoir_properties = createReservoirProperties(
        //! \ogs_file_param{prj__processes__process__WELLBORE_SIMULATOR__reservoir_properties}
        config.getConfigSubtree("reservoir_properties"), parameters);

    auto const& productivity_index = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__WELLBORE_SIMULATOR__productivity_index}
        "productivity_index", parameters, 1, &mesh);
    DBUG("Use '{:s}' as productivity_index parameter.",
         productivity_index.name);

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    DBUG("Check the media properties of WellboreSimulator process ...");
    checkMPLProperties(mesh, media_map);
    DBUG("Media properties verified.");

    WellboreSimulatorProcessData process_data{std::move(media_map),
                                              specific_body_force,
                                              std::move(wellbore_geometry),
                                              well_ref_pressure,
                                              well_ref_enthalpy,
                                              std::move(reservoir_properties),
                                              productivity_index,
                                              heat_exchange_with_formation,
                                              has_gravity};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<WellboreSimulatorProcess>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables));
}

}  // namespace WellboreSimulator
}  // namespace ProcessLib
