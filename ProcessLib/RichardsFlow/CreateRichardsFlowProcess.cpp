/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateRichardsFlowProcess.h"

#include "ProcessLib/Output/ParseSecondaryVariables.h"
#include "ProcessLib/Parameter/ConstantParameter.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "CreateRichardsFlowMaterialProperties.h"
#include "RichardsFlowProcess.h"
#include "RichardsFlowProcessData.h"

namespace ProcessLib
{
namespace RichardsFlow
{
std::unique_ptr<Process> createRichardsFlowProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "RICHARDS_FLOW");

    DBUG("Create RichardsFlowProcess.");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__RICHARDS_FLOW__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    auto per_process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__RICHARDS_FLOW__process_variables__process_variable}
         "process_variable"});
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    process_variables.push_back(std::move(per_process_variables));

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"RichardsFlow_pressure"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    // Specific body force
    std::vector<double> const b =
        //! \ogs_file_param{prj__processes__process__RICHARDS_FLOW__specific_body_force}
        config.getConfigParameter<std::vector<double>>("specific_body_force");
    assert(b.size() > 0 && b.size() < 4);
    Eigen::VectorXd specific_body_force(b.size());
    bool const has_gravity = MathLib::toVector(b).norm() > 0;
    if (has_gravity)
        std::copy_n(b.data(), b.size(), specific_body_force.data());

    // has mass lumping
    //! \ogs_file_param{prj__processes__process__RICHARDS_FLOW__mass_lumping}
    auto mass_lumping = config.getConfigParameter<bool>("mass_lumping");

    auto& temperature = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__RICHARDS_FLOW__temperature}
        "temperature", parameters, 1);

    //! \ogs_file_param{prj__processes__process__RICHARDS_FLOW__material_property}
    auto const& mat_config = config.getConfigSubtree("material_property");

    boost::optional<MeshLib::PropertyVector<int> const&> material_ids;
    if (mesh.getProperties().existsPropertyVector<int>("MaterialIDs"))
    {
        INFO("The Richards flow is in heterogeneous porous media.");
        material_ids =
            *mesh.getProperties().getPropertyVector<int>("MaterialIDs");
    }
    else
    {
        INFO("The Richards flow is in homogeneous porous media.");
    }
    std::unique_ptr<RichardsFlowMaterialProperties> material =
        createRichardsFlowMaterialProperties(mat_config, material_ids,
                                             parameters);
    RichardsFlowProcessData process_data{std::move(material),
                                         specific_body_force, has_gravity,
                                         mass_lumping, temperature};

    return std::make_unique<RichardsFlowProcess>(
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller),
        mat_config, curves);
}

}  // namespace RichardsFlow
}  // namespace ProcessLib
