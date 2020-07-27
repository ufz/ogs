/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateSmallDeformationNonlocalProcess.h"

#include <cassert>

#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "SmallDeformationNonlocalProcess.h"
#include "SmallDeformationNonlocalProcessData.h"

namespace ProcessLib
{
namespace SmallDeformationNonlocal
{
template <int DisplacementDim>
std::unique_ptr<Process> createSmallDeformationNonlocalProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "SMALL_DEFORMATION_NONLOCAL");
    DBUG("Create SmallDeformationNonlocalProcess.");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION_NONLOCAL__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    auto per_process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__SMALL_DEFORMATION_NONLOCAL__process_variables__process_variable}
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

    auto solid_constitutive_relations =
        MaterialLib::Solids::createConstitutiveRelations<DisplacementDim>(
            parameters, local_coordinate_system, config);

    // Solid density
    auto& solid_density = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__SMALL_DEFORMATION_NONLOCAL__solid_density}
        "solid_density", parameters, 1, &mesh);
    DBUG("Use '{:s}' as solid density parameter.", solid_density.name);

    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION_NONLOCAL__specific_body_force}
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

    // Reference temperature
    const auto& reference_temperature =
        //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION_NONLOCAL__reference_temperature}
        config.getConfigParameter<double>(
            "reference_temperature", std::numeric_limits<double>::quiet_NaN());

    auto const internal_length =
        //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION_NONLOCAL__internal_length}
        config.getConfigParameter<double>("internal_length");

    SmallDeformationNonlocalProcessData<DisplacementDim> process_data{
        materialIDs(mesh),     std::move(solid_constitutive_relations),
        solid_density,         specific_body_force,
        reference_temperature, internal_length * internal_length};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<SmallDeformationNonlocalProcess<DisplacementDim>>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables));
}

template std::unique_ptr<Process> createSmallDeformationNonlocalProcess<2>(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<Process> createSmallDeformationNonlocalProcess<3>(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib
