// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateHTProcess.h"

#include "HTLocalAssemblerInterface.h"
#include "HTProcess.h"
#include "HTProcessData.h"
#include "MaterialLib/MPL/CheckMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Utils/GetElementRotationMatrices.h"
#include "MeshLib/Utils/GetSpaceDimension.h"
#include "NumLib/NumericalStability/CreateNumericalStabilization.h"
#include "ParameterLib/ConstantParameter.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Common/HydraulicProcess/checkVolumeBalanceEquationSetting.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/SurfaceFlux/SurfaceFluxData.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
namespace HT
{

void checkThermalExpansivitySetting(
    MaterialPropertyLib::MaterialSpatialDistributionMap const& media_map)
{
    for (auto const& medium : media_map.media())
    {
        auto const& solid_phase =
            medium->phase(MaterialPropertyLib::PhaseName::Solid);

        bool const has_thermal_expansivity = solid_phase.hasProperty(
            MaterialPropertyLib::PropertyType::thermal_expansivity);
        if (has_thermal_expansivity)
        {
            bool const has_biot_constant = medium->hasProperty(
                MaterialPropertyLib::PropertyType::biot_coefficient);
            if (!has_biot_constant)
            {
                OGS_FATAL(
                    "Since the solid phase has thermal expansivity, it must "
                    "also have the Biot constant. Please add the property "
                    "'biot_coefficient' to the `properties` in the material "
                    "configuration.");
            }
        }
    }
}
void checkMPLProperties(
    MeshLib::Mesh const& mesh,
    MaterialPropertyLib::MaterialSpatialDistributionMap const& media_map)
{
    std::array const required_property_medium = {
        MaterialPropertyLib::PropertyType::permeability,
        MaterialPropertyLib::PropertyType::porosity,
        MaterialPropertyLib::PropertyType::thermal_conductivity,
        MaterialPropertyLib::PropertyType::thermal_longitudinal_dispersivity,
        MaterialPropertyLib::PropertyType::thermal_transversal_dispersivity};

    std::array const required_property_liquid_phase = {
        MaterialPropertyLib::PropertyType::viscosity,
        MaterialPropertyLib::PropertyType::density,
        MaterialPropertyLib::PropertyType::specific_heat_capacity,
        MaterialPropertyLib::PropertyType::thermal_conductivity};

    std::array const required_property_solid_phase = {
        MaterialPropertyLib::PropertyType::specific_heat_capacity,
        MaterialPropertyLib::PropertyType::density,
        MaterialPropertyLib::PropertyType::thermal_conductivity,
        MaterialPropertyLib::PropertyType::storage};

    std::array<MaterialPropertyLib::PropertyType, 0> const
        required_gas_properties{};

    MaterialPropertyLib::checkMaterialSpatialDistributionMap(
        mesh, media_map, required_property_medium,
        required_property_solid_phase, required_property_liquid_phase,
        required_gas_properties);

    checkThermalExpansivitySetting(media_map);
}

std::unique_ptr<Process> createHTProcess(
    std::string const& name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "HT");

    DBUG("Create HTProcess.");

    auto const coupling_scheme =
        //! \ogs_file_param{prj__processes__process__HT__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(coupling_scheme && (*coupling_scheme == "staggered"));

    /// \section processvariablesht Process Variables

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
        /// Primary process variables as they appear in the global component
        /// vector:
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

    /// \section parametersht Process Parameters
    std::vector<double> const b =
        //! \ogs_file_param{prj__processes__process__HT__specific_body_force}
        config.getConfigParameter<std::vector<double>>("specific_body_force");
    assert(!b.empty() && b.size() < 4);
    int const mesh_space_dimension =
        MeshLib::getSpaceDimension(mesh.getNodes());
    if (static_cast<int>(b.size()) != mesh_space_dimension)
    {
        OGS_FATAL(
            "specific body force (gravity vector) has {:d} components, "
            "mesh dimension is {:d}",
            b.size(), mesh_space_dimension);
    }

    // Specific body force parameter.
    Eigen::VectorXd specific_body_force(b.size());
    bool const has_gravity = MathLib::toVector(b).norm() > 0;
    if (has_gravity)
    {
        std::copy_n(b.data(), b.size(), specific_body_force.data());
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
    checkMPLProperties(mesh, media_map);
    DBUG("Media properties verified.");

    auto stabilizer = NumLib::createNumericalStabilization(mesh, config);

    auto const* aperture_size_parameter = &ParameterLib::findParameter<double>(
        ProcessLib::Process::constant_one_parameter_name, parameters, 1);
    auto const aperture_config =
        //! \ogs_file_param{prj__processes__process__HT__aperture_size}
        config.getConfigSubtreeOptional("aperture_size");
    if (aperture_config)
    {
        aperture_size_parameter = &ParameterLib::findParameter<double>(
            //! \ogs_file_param_special{prj__processes__process__HT__aperture_size__parameter}
            *aperture_config, "parameter", parameters, 1);
    }

    auto const rotation_matrices = MeshLib::getElementRotationMatrices(
        mesh_space_dimension, mesh.getDimension(), mesh.getElements());
    std::vector<Eigen::VectorXd> projected_specific_body_force_vectors;
    projected_specific_body_force_vectors.reserve(rotation_matrices.size());

    std::transform(rotation_matrices.begin(), rotation_matrices.end(),
                   std::back_inserter(projected_specific_body_force_vectors),
                   [&specific_body_force](const auto& R)
                   { return R * R.transpose() * specific_body_force; });

    auto const equation_balance_type_str =
        //! \ogs_file_param{prj__processes__process__HT__equation_balance_type}
        config.getConfigParameter<std::string>("equation_balance_type",
                                               "volume");
    if (equation_balance_type_str != "volume" &&
        equation_balance_type_str != "mass")
    {
        OGS_FATAL(
            "Invalid equation_balance_type '{}'. Supported values: 'volume' or "
            "'mass'.",
            equation_balance_type_str);
    }
    bool const is_volume_balance_equation_type =
        (equation_balance_type_str == "volume");

    if (is_volume_balance_equation_type)
    {
        ProcessLib::Common::HydraulicProcess::checkVolumeBalanceEquationSetting(
            media_map);
    }

    HTProcessData process_data{std::move(media_map),
                               has_gravity,
                               heat_transport_process_id,
                               hydraulic_process_id,
                               std::move(stabilizer),
                               projected_specific_body_force_vectors,
                               mesh_space_dimension,
                               *aperture_size_parameter,
                               is_volume_balance_equation_type,
                               NumLib::ShapeMatrixCache{integration_order}};

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
