/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "UniformRobinBoundaryCondition.h"
#include "BoundaryConditionConfig.h"
#include "MeshGeoToolsLib/BoundaryElementsSearcher.h"

namespace ProcessLib
{
std::unique_ptr<UniformRobinBoundaryCondition>
createUniformRobinBoundaryCondition(
    BoundaryConditionConfig const& config,
    MeshGeoToolsLib::BoundaryElementsSearcher& searcher,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    unsigned const integration_order, const unsigned global_dim)
{
    DBUG("Constructing RobinBcConfig from config.");
    //! \ogs_file_param{boundary_condition__type}
    config.config.checkConfigParameter("type", "UniformRobin");

    //! \ogs_file_param{boundary_condition__UniformRobin__alpha}
    double const alpha = config.config.getConfigParameter<double>("alpha");
    //! \ogs_file_param{boundary_condition__UniformRobin__f_0}
    double const f_0 = config.config.getConfigParameter<double>("f_0");

    auto& elems = searcher.getBoundaryElements(config.geometry);

    return std::unique_ptr<UniformRobinBoundaryCondition>(
        new UniformRobinBoundaryCondition(
            integration_order, dof_table, variable_id, config.component_id,
            global_dim, elems, UniformRobinBoundaryConditionData{alpha, f_0}));
}

}  // ProcessLib
