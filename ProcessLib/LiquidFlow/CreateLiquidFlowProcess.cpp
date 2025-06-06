/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * Created on August 19, 2016, 1:30 PM
 */
#include "CreateLiquidFlowProcess.h"

#include <algorithm>
#include <typeinfo>

#include "LiquidFlowProcess.h"
#include "MaterialLib/MPL/CheckMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Properties/Constant.h"
#include "MaterialLib/MPL/Utils/CheckMPLPhasesForSinglePhaseFlow.h"
#include "MeshLib/Utils/GetElementRotationMatrices.h"
#include "MeshLib/Utils/GetSpaceDimension.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
namespace LiquidFlow
{
void checkMPLProperties(
    MeshLib::Mesh const& mesh,
    MaterialPropertyLib::MaterialSpatialDistributionMap const& media_map,
    bool const is_equation_type_volume)
{
    MaterialPropertyLib::checkMPLPhasesForSinglePhaseFlow(mesh, media_map);

    std::array const required_medium_properties = {
        MaterialPropertyLib::reference_temperature,
        MaterialPropertyLib::PropertyType::permeability,
        MaterialPropertyLib::PropertyType::porosity,
        MaterialPropertyLib::PropertyType::storage};

    std::array const required_liquid_properties = {
        MaterialPropertyLib::PropertyType::viscosity,
        MaterialPropertyLib::PropertyType::density};

    // Check Constant-type density.
    if (is_equation_type_volume)
    {
        for (auto const& medium : media_map.media())
        {
            // auto const& medium = *media_map.getMedium(element_id);
            auto const& fluid_phase_density =
                fluidPhase(*medium)[MaterialPropertyLib::PropertyType::density];
            if (typeid(fluid_phase_density) !=
                typeid(MaterialPropertyLib::Constant))
            {
                OGS_FATAL(
                    "Since `equation_balance_type` is set to `volume`,the "
                    "phase density type must be `Constant`. Note: by "
                    "default, `equation_balance_type` is set to `volume`.");
            }
        }
    }

    for (auto const& medium : media_map.media())
    {
        checkRequiredProperties(*medium, required_medium_properties);
        checkRequiredProperties(fluidPhase(*medium),
                                required_liquid_properties);
    }
    DBUG("Media properties verified.");
}

EquationBalanceType covertEquationBalanceTypeFromString(
    std::string_view const type_in_string)
{
    if (type_in_string == "volume")
    {
        return EquationBalanceType::volume;
    }
    if (type_in_string == "mass")
    {
        return EquationBalanceType::mass;
    }
    OGS_FATAL(
        "The value of `equation_balance_type` must be either `volume` or "
        "`mass`");
};

std::unique_ptr<Process> createLiquidFlowProcess(
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
    config.checkConfigParameter("type", "LIQUID_FLOW");

    DBUG("Create LiquidFlowProcess.");

    /// \section processvariableslf Process Variables

    //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    /// Primary process variables as they appear in the global component vector:
    auto per_process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__LIQUID_FLOW__process_variables__process_variable}
         "process_variable"});
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    process_variables.push_back(std::move(per_process_variables));

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    /// \section parameterslf Process Parameters
    std::vector<double> const b =
        //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__specific_body_force}
        config.getConfigParameter<std::vector<double>>("specific_body_force");
    int const mesh_space_dimension =
        MeshLib::getSpaceDimension(mesh.getNodes());
    if (static_cast<int>(b.size()) != mesh_space_dimension)
    {
        OGS_FATAL(
            "specific body force (gravity vector) has {:d} components, but the "
            "space dimension is {:d}.",
            b.size(), mesh_space_dimension);
    }

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

    auto const is_linear =
        //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__linear}
        config.getConfigParameter<bool>("linear", false);
    if (is_linear)
    {
        INFO("LiquidFlow process is set to be linear.");
    }

    auto const equation_balance_type_str =
        //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__equation_balance_type}
        config.getConfigParameter<std::string>("equation_balance_type",
                                               "volume");

    DBUG("Check the media properties of LiquidFlow process ...");
    checkMPLProperties(mesh, media_map, equation_balance_type_str == "volume");
    DBUG("Media properties verified.");

    auto const* aperture_size_parameter = &ParameterLib::findParameter<double>(
        ProcessLib::Process::constant_one_parameter_name, parameters, 1);

    auto const aperture_config =
        //! \ogs_file_param{prj__processes__process__LIQUID_FLOW__aperture_size}
        config.getConfigSubtreeOptional("aperture_size");
    if (aperture_config)
    {
        aperture_size_parameter = &ParameterLib::findParameter<double>(
            //! \ogs_file_param_special{prj__processes__process__LIQUID_FLOW__aperture_size__parameter}
            *aperture_config, "parameter", parameters, 1);
    }

    // The uniqueness of phase has already been checked in
    // `checkMPLProperties`.
    MaterialPropertyLib::Variable const phase_variable =
        (*ranges::begin(media_map.media()))->hasPhase("Gas")
            ? MaterialPropertyLib::Variable::gas_phase_pressure
            : MaterialPropertyLib::Variable::liquid_phase_pressure;

    LiquidFlowData process_data{
        covertEquationBalanceTypeFromString(equation_balance_type_str),
        std::move(media_map),
        MeshLib::getElementRotationMatrices(
            mesh_space_dimension, mesh.getDimension(), mesh.getElements()),
        mesh_space_dimension,
        std::move(specific_body_force),
        has_gravity,
        *aperture_size_parameter,
        NumLib::ShapeMatrixCache{integration_order},
        phase_variable};

    return std::make_unique<LiquidFlowProcess>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        std::move(surfaceflux), is_linear);
}

}  // namespace LiquidFlow
}  // namespace ProcessLib
