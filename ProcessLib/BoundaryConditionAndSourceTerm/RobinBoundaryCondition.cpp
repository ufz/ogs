/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "RobinBoundaryCondition.h"

#include "ParameterLib/Utils.h"

namespace ProcessLib
{
RobinBoundaryConditionConfig parseRobinBoundaryCondition(
    BaseLib::ConfigTree const& config)
{
    DBUG("parsing RobinBoundaryConditionConfig.");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "Robin");

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__Robin__alpha}
    auto const alpha_name = config.getConfigParameter<std::string>("alpha");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__Robin__u_0}
    auto const u_0_name = config.getConfigParameter<std::string>("u_0");

    auto const area_parameter_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__Robin__area_parameter}
        config.getConfigParameterOptional<std::string>("area_parameter");
    if (area_parameter_name.has_value())
    {
        DBUG("area parameter name '{:s}'", __FUNCTION__,
             area_parameter_name.value());
    }

    return {alpha_name, u_0_name, area_parameter_name};
}

std::unique_ptr<RobinBoundaryCondition> createRobinBoundaryCondition(
    RobinBoundaryConditionConfig const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const shapefunction_order, unsigned const global_dim,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    DBUG("Constructing RobinBcConfig.");

    if (bc_mesh.getDimension() >= global_dim)
    {
        OGS_FATAL(
            "The dimension ({:d}) of the given boundary mesh '{:s}' is not "
            "lower than the bulk dimension ({:d}).",
            bc_mesh.getDimension(), bc_mesh.getName(), global_dim);
    }

    auto const& alpha = ParameterLib::findParameter<double>(
        config.alpha_name, parameters, 1, &bc_mesh);
    auto const& u_0 = ParameterLib::findParameter<double>(
        config.u_0_name, parameters, 1, &bc_mesh);

    ParameterLib::Parameter<double> const* integral_measure(nullptr);
    if (global_dim - bc_mesh.getDimension() != 1)
    {
        if (!config.area_parameter_name)
        {
            OGS_FATAL("{}: tag <area_parameter> required in Robin BC.",
                      __FUNCTION__);
        }
        integral_measure = &ParameterLib::findParameter<double>(
            config.area_parameter_name.value(), parameters, 1, &bc_mesh);
    }

    // In case of partitioned mesh the boundary could be empty, i.e. there is no
    // boundary condition.
#ifdef USE_PETSC
    // This can be extracted to createBoundaryCondition() but then the config
    // parameters are not read and will cause an error.
    // TODO (naumov): Add a function to ConfigTree for skipping the tags of the
    // subtree and move the code up in createBoundaryCondition().
    if (bc_mesh.getDimension() == 0 && bc_mesh.getNumberOfNodes() == 0 &&
        bc_mesh.getNumberOfElements() == 0)
    {
        return nullptr;
    }
#endif  // USE_PETSC

    return std::make_unique<RobinBoundaryCondition>(
        integration_order, shapefunction_order, dof_table, variable_id,
        component_id, global_dim, bc_mesh,
        RobinBoundaryConditionData{alpha, u_0, integral_measure});
}

}  // namespace ProcessLib
