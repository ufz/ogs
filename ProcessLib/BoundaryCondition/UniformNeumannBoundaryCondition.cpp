/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "UniformNeumannBoundaryCondition.h"
#include "MeshGeoToolsLib/BoundaryElementsSearcher.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "BoundaryConditionConfig.h"

namespace ProcessLib
{
std::unique_ptr<UniformNeumannBoundaryCondition>
createUniformNeumannBoundaryCondition(
    BoundaryConditionConfig const& config,
    MeshGeoToolsLib::BoundaryElementsSearcher& searcher,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    unsigned const integration_order, const unsigned global_dim)
{
    DBUG("Constructing NeumannBcConfig from config.");
    //! \ogs_file_param{boundary_condition__type}
    config.config.checkConfigParameter("type", "UniformNeumann");

    //! \ogs_file_param{boundary_condition__UniformNeumann__value}
    double const value = config.config.getConfigParameter<double>("value");
    DBUG("Using value %g", value);

    auto& elems = searcher.getBoundaryElements(config.geometry);

    return std::unique_ptr<UniformNeumannBoundaryCondition>(
        new UniformNeumannBoundaryCondition(integration_order, dof_table,
                                            variable_id, config.component_id,
                                            global_dim, elems, value));
}

}  // ProcessLib
