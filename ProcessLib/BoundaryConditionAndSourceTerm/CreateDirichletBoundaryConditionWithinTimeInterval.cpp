/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BaseLib/ConfigTree.h"
#include "BaseLib/TimeInterval.h"
#include "DirichletBoundaryConditionWithinTimeInterval.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/Utils.h"

namespace ProcessLib
{
std::unique_ptr<BoundaryCondition>
createDirichletBoundaryConditionWithinTimeInterval(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    int const component_id,
    const std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters)
{
    DBUG(
        "Constructing DirichletBoundaryConditionWithinTimeInterval from "
        "config.");

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "DirichletWithinTimeInterval");

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__DirichletWithinTimeInterval__parameter}
    auto const param_name = config.getConfigParameter<std::string>("parameter");
    DBUG("Using parameter {:s}", param_name);

    auto& param = ParameterLib::findParameter<double>(param_name, parameters, 1,
                                                      &bc_mesh);

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__DirichletWithinTimeInterval__time_interval}
    config.peekConfigParameter<std::string>("time_interval");
    auto time_interval = BaseLib::createTimeInterval(config);

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

    return std::make_unique<DirichletBoundaryConditionWithinTimeInterval>(
        std::move(time_interval), param, bc_mesh, dof_table_bulk, variable_id,
        component_id);
}
}  // namespace ProcessLib
