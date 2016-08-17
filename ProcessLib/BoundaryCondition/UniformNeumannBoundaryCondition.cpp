/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "UniformNeumannBoundaryCondition.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
std::unique_ptr<UniformNeumannBoundaryCondition>
createUniformNeumannBoundaryCondition(
    BaseLib::ConfigTree const& config,
    std::vector<MeshLib::Element*>&& elements,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const global_dim,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters)
{
    DBUG("Constructing NeumannBcConfig from config.");
    //! \ogs_file_param{boundary_condition__type}
    config.checkConfigParameter("type", "UniformNeumann");

    //! \ogs_file_param{boundary_condition__UniformNeumann__parameter}
    auto const param_name = config.getConfigParameter<std::string>("parameter");
    DBUG("Using parameter %g", param_name);

    auto& param = findParameter<double>(param_name, parameters, 1);

    return std::unique_ptr<UniformNeumannBoundaryCondition>(
        new UniformNeumannBoundaryCondition(
            integration_order, dof_table, variable_id, component_id,
            global_dim, std::move(elements), param));
}

}  // ProcessLib
