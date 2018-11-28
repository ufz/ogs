/**
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * File:   DirichletBoundaryConditionAuxiliaryFunctions.cpp
 *
 * Created on November 28, 2018, 11:26 AM
 */
#include "DirichletBoundaryConditionAuxiliaryFunctions.h"

#include "NumLib/IndexValueVector.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/Parameter/Parameter.h"

namespace ProcessLib
{
void checkParametersOfDirichletBoundaryCondition(
    MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
    int const variable_id,
    int const component_id)
{
    if (variable_id >=
            static_cast<int>(dof_table_bulk.getNumberOfVariables()) ||
        component_id >=
            dof_table_bulk.getNumberOfVariableComponents(variable_id))
    {
        OGS_FATAL(
            "Variable id or component id too high. Actual values: (%d, "
            "%d), maximum values: (%d, %d).",
            variable_id, component_id, dof_table_bulk.getNumberOfVariables(),
            dof_table_bulk.getNumberOfVariableComponents(variable_id));
    }

    if (!bc_mesh.getProperties().existsPropertyVector<std::size_t>(
            "bulk_node_ids"))
    {
        OGS_FATAL(
            "The required bulk node ids map does not exist in the boundary "
            "mesh '%s'.",
            bc_mesh.getName().c_str());
    }

    DBUG(
        "Found %d nodes for Dirichlet BCs for the variable %d and "
        "component "
        "%d",
        bc_mesh.getNodes().size(), variable_id, component_id);
}

void getEssentialBCValuesLocal(
    Parameter<double> const& parameter, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
    int const variable_id, int const component_id, const double t,
    GlobalVector const& /*x*/,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values)
{
    SpatialPosition pos;

    bc_values.ids.clear();
    bc_values.values.clear();

    // convert mesh node ids to global index for the given component
    bc_values.ids.reserve(bc_values.ids.size() + bc_mesh.getNumberOfNodes());
    bc_values.values.reserve(bc_values.values.size() +
                             bc_mesh.getNumberOfNodes());
    for (auto const* const node : bc_mesh.getNodes())
    {
        auto const id = node->getID();
        pos.setNodeID(node->getID());
        // TODO: that might be slow, but only done once
        auto const global_index = dof_table_boundary.getGlobalIndex(
            {bc_mesh.getID(), MeshLib::MeshItemType::Node, id}, variable_id,
            component_id);
        if (global_index == NumLib::MeshComponentMap::nop)
            continue;
        // For the DDC approach (e.g. with PETSc option), the negative
        // index of global_index means that the entry by that index is a ghost
        // one, which should be dropped. Especially for PETSc routines
        // MatZeroRows and MatZeroRowsColumns, which are called to apply the
        // Dirichlet BC, the negative index is not accepted like other matrix or
        // vector PETSc routines. Therefore, the following if-condition is
        // applied.
        if (global_index >= 0)
        {
            bc_values.ids.emplace_back(global_index);
            bc_values.values.emplace_back(parameter(t, pos).front());
        }
    }
}
} // end of name space
