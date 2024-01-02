/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateLargeDeformationProcess.h"

#include <cassert>

#include "ConstitutiveRelations/CreateConstitutiveSetting.h"
#include "LargeDeformationProcess.h"
#include "LargeDeformationProcessData.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
namespace LargeDeformation
{
template <int DisplacementDim>
std::unique_ptr<Process> createLargeDeformationProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "LARGE_DEFORMATION");
    DBUG("Create LargeDeformationProcess.");

    /// \section processvariablesld Process Variables

    //! \ogs_file_param{prj__processes__process__LARGE_DEFORMATION__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    /// Primary process variables as they appear in the global component vector:
    auto per_process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__LARGE_DEFORMATION__process_variables__process_variable}
         "process_variable"});

    DBUG("Associate displacement with process variable '{:s}'.",
         per_process_variables.back().get().getName());

    if (per_process_variables.back().get().getNumberOfGlobalComponents() !=
        DisplacementDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '{:s}' is different "
            "from the displacement dimension: got {:d}, expected {:d}",
            per_process_variables.back().get().getName(),
            per_process_variables.back().get().getNumberOfGlobalComponents(),
            DisplacementDim);
    }
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    process_variables.push_back(std::move(per_process_variables));

    /// \section parametersld Process Parameters
    auto solid_constitutive_relations =
        ConstitutiveRelations::CreateConstitutiveSetting<DisplacementDim>::
            createSolidConstitutiveRelations(parameters,
                                             local_coordinate_system, config);

    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__LARGE_DEFORMATION__specific_body_force}
            config.getConfigParameter<std::vector<double>>(
                "specific_body_force");
        if (b.size() != DisplacementDim)
        {
            OGS_FATAL(
                "The size of the specific body force vector does not match the "
                "displacement dimension. Vector size is {:d}, displacement "
                "dimension is {:d}",
                b.size(), DisplacementDim);
        }

        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    // Reference temperature
    auto const reference_temperature = ParameterLib::findOptionalTagParameter<
        double>(
        //! \ogs_file_param_special{prj__processes__process__LARGE_DEFORMATION__reference_temperature}
        config, "reference_temperature", parameters, 1, &mesh);
    if (reference_temperature)
    {
        DBUG("Use '{:s}' as reference temperature parameter.",
             (*reference_temperature).name);
    }

    // Initial stress conditions
    auto const initial_stress = ParameterLib::findOptionalTagParameter<double>(
        //! \ogs_file_param_special{prj__processes__process__LARGE_DEFORMATION__initial_stress}
        config, "initial_stress", parameters,
        // Symmetric tensor size, 4 or 6, not a Kelvin vector.
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim),
        &mesh);

    LargeDeformationProcessData<DisplacementDim> process_data{
        materialIDs(mesh),
        std::move(media_map),
        std::move(solid_constitutive_relations),
        initial_stress,
        specific_body_force,
        reference_temperature};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<LargeDeformationProcess<DisplacementDim>>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables));
}

template std::unique_ptr<Process> createLargeDeformationProcess<2>(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media);

template std::unique_ptr<Process> createLargeDeformationProcess<3>(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media);

}  // namespace LargeDeformation
}  // namespace ProcessLib
