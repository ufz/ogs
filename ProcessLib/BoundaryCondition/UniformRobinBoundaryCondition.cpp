/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "UniformRobinBoundaryCondition.h"

namespace ProcessLib
{
std::unique_ptr<UniformRobinBoundaryCondition>
createUniformRobinBoundaryCondition(
    BaseLib::ConfigTree const& config,
    std::vector<MeshLib::Element*>&& elements,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const global_dim)
{
    DBUG("Constructing RobinBcConfig from config.");
    //! \ogs_file_param{boundary_condition__type}
    config.checkConfigParameter("type", "UniformRobin");

    //! \ogs_file_param{boundary_condition__UniformRobin__alpha}
    double const alpha = config.getConfigParameter<double>("alpha");
    //! \ogs_file_param{boundary_condition__UniformRobin__f_0}
    double const f_0 = config.getConfigParameter<double>("f_0");

    return std::unique_ptr<UniformRobinBoundaryCondition>(
        new UniformRobinBoundaryCondition(
            integration_order, dof_table, variable_id, component_id, global_dim,
            std::move(elements),
            UniformRobinBoundaryConditionData{alpha, f_0}));
}

}  // ProcessLib
