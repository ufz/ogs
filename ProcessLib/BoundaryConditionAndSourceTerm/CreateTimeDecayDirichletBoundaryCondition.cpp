/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on 2025-07-29 14:36:50
 */

#include "CreateTimeDecayDirichletBoundaryCondition.h"

#include <numeric>

#include "BaseLib/ConfigTree.h"
#include "BoundaryCondition.h"
#include "MeshLib/Mesh.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/Utils.h"
#include "TimeDecayDirichletBoundaryCondition.h"

namespace ProcessLib
{
std::unique_ptr<BoundaryCondition> createTimeDecayDirichletBoundaryCondition(
    int const variable_id, int const component_id,
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    DBUG("Create TimeDecayDirichlet.");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "TimeDecayDirichlet");

    auto const parameter_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__TimeDecayDirichlet__time_decay_parameter}
        config.getConfigParameter<std::string>("time_decay_parameter");
    DBUG("Using parameter {:s}", parameter_name);

    auto const& time_decay_parameter = ParameterLib::findParameter<double>(
        parameter_name, parameters, 1, &bc_mesh);

    auto const lower_limit =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__TimeDecayDirichlet__lower_limit}
        config.getConfigParameter<double>("lower_limit");

    // In case of partitioned mesh the boundary could be empty, i.e. there
    // is no boundary condition.
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

    return std::make_unique<TimeDecayDirichletBoundaryCondition>(
        variable_id, component_id, bc_mesh, dof_table_bulk,
        time_decay_parameter, lower_limit);
}

}  // namespace ProcessLib
