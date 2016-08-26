/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NeumannBoundaryCondition.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
std::unique_ptr<NeumannBoundaryCondition>
createNeumannBoundaryCondition(
    BaseLib::ConfigTree const& config,
    std::vector<MeshLib::Element*>&& elements,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const global_dim,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters)
{
    DBUG("Constructing Neumann BC from config.");
    //! \ogs_file_param{boundary_condition__type}
    config.checkConfigParameter("type", "Neumann");

    //! \ogs_file_param{boundary_condition__Neumann__parameter}
    auto const param_name = config.getConfigParameter<std::string>("parameter");
    DBUG("Using parameter %s", param_name.c_str());

    auto const& param = findParameter<double>(param_name, parameters, 1);

    return std::unique_ptr<NeumannBoundaryCondition>(
        new NeumannBoundaryCondition(
            integration_order, dof_table, variable_id, component_id,
            global_dim, std::move(elements), param));
}

}  // ProcessLib
