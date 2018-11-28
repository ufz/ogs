/**
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * File:   DirichletBoundaryConditionWithinTimeInterval.cpp
 *
 * Created on November 26, 2018, 4:59 PM
 */
#include "DirichletBoundaryConditionWithinTimeInterval.h"

#include "BaseLib/ConfigTree.h"

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"
#include "NumLib/TimeStepping/TimeInterval.h"

#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
DirichletBoundaryConditionWithinTimeInterval::
    DirichletBoundaryConditionWithinTimeInterval(
        std::unique_ptr<NumLib::TimeInterval> time_interval,
        Parameter<double> const& parameter, MeshLib::Mesh const& bc_mesh,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id)
    : DirichletBoundaryCondition(parameter, bc_mesh, dof_table_bulk,
                                 variable_id, component_id),
      _time_interval(std::move(time_interval))
{
}

void DirichletBoundaryConditionWithinTimeInterval::getEssentialBCValues(
    const double t, GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    if (_time_interval->contains(t))
    {
        getEssentialBCValuesLocal(t, x, bc_values);
        return;
    }

    bc_values.ids.clear();
    bc_values.values.clear();
}

std::unique_ptr<DirichletBoundaryCondition>
createDirichletBoundaryConditionWithinTimeInterval(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    int const component_id,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters)
{
    DBUG(
        "Constructing DirichletBoundaryConditionWithinTimeInterval from "
        "config.");

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "DirichletWithinTimeInterval");

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__DirichletWithinTimeInterval__parameter}
    auto const param_name = config.getConfigParameter<std::string>("parameter");
    DBUG("Using parameter %s", param_name.c_str());

    auto& param = findParameter<double>(param_name, parameters, 1);

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

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__DirichletWithinTimeInterval__time_interval}
    config.peekConfigParameter<std::string>("time_interval");

    return std::make_unique<DirichletBoundaryConditionWithinTimeInterval>(
        NumLib::createTimeInterval(config), param, bc_mesh,
        dof_table_bulk, variable_id, component_id);
}

}  // namespace ProcessLib
