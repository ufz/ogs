/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * File:   DirichletBoundaryConditionAuxiliaryFunctions.cpp
 *
 * Created on November 28, 2018, 11:26 AM
 */
#include "DirichletBoundaryConditionAuxiliaryFunctions.h"

#include <range/v3/range/conversion.hpp>

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"
#include "ParameterLib/Parameter.h"

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
            "Variable id or component id too high. Actual values: ({:d}, "
            "{:d}), maximum values: ({:d}, {:d}).",
            variable_id, component_id, dof_table_bulk.getNumberOfVariables(),
            dof_table_bulk.getNumberOfVariableComponents(variable_id));
    }

    if (!bc_mesh.getProperties().hasPropertyVector(
            std::string(MeshLib::getBulkIDString(MeshLib::MeshItemType::Node))))
    {
        OGS_FATAL(
            "The required bulk node ids map does not exist in the boundary "
            "mesh '{:s}'.",
            bc_mesh.getName());
    }

    if (!bc_mesh.getProperties().existsPropertyVector<std::size_t>(
            MeshLib::getBulkIDString(MeshLib::MeshItemType::Node)))
    {
        OGS_FATAL(
            "The required bulk node ids map exists in the boundary mesh '{:s}' "
            "but has wrong data type (should be equivalent to C++ data type "
            "std::size_t which is an unsigned integer of size {:d} or UInt64 "
            "in vtk terminology).",
            bc_mesh.getName(), sizeof(std::size_t));
    }

    DBUG(
        "Found {:d} nodes for Dirichlet BCs for the variable {:d} and "
        "component {:d}",
        bc_mesh.getNodes().size(), variable_id, component_id);
}

void getEssentialBCValuesLocal(
    ParameterLib::Parameter<double> const& parameter,
    MeshLib::Mesh const& bc_mesh,
    std::vector<std::size_t> const& nodes_in_bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
    int const variable_id, int const component_id, const double t,
    GlobalVector const& /*x*/,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values)
{
    ParameterLib::SpatialPosition pos;

    bc_values.ids.clear();
    bc_values.values.clear();

    // convert mesh node ids to global index for the given component
    bc_values.ids.reserve(nodes_in_bc_mesh.size());
    bc_values.values.reserve(nodes_in_bc_mesh.size());
    for (auto const node_id : nodes_in_bc_mesh)
    {
        // TODO: that might be slow, but only done once
        auto const global_index = dof_table_boundary.getGlobalIndex(
            {bc_mesh.getID(), MeshLib::MeshItemType::Node, node_id},
            variable_id, component_id);
        if (global_index == NumLib::MeshComponentMap::nop)
        {
            continue;
        }
        // For the DDC approach (e.g. with PETSc option), the negative
        // index of global_index means that the entry by that index is a ghost
        // one, which should be dropped. Especially for PETSc routines
        // MatZeroRows and MatZeroRowsColumns, which are called to apply the
        // Dirichlet BC, the negative index is not accepted like other matrix or
        // vector PETSc routines. Therefore, the following if-condition is
        // applied.
        if (global_index >= 0)
        {
            pos.setNodeID(node_id);
            pos.setCoordinates(*bc_mesh.getNode(node_id));
            bc_values.ids.emplace_back(global_index);
            bc_values.values.emplace_back(parameter(t, pos).front());
        }
    }
}

void getEssentialBCValuesLocal(
    ParameterLib::Parameter<double> const& parameter,
    MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
    int const variable_id, int const component_id, const double t,
    GlobalVector const& x, NumLib::IndexValueVector<GlobalIndexType>& bc_values)
{
    getEssentialBCValuesLocal(
        parameter, bc_mesh,
        bc_mesh.getNodes() | MeshLib::views::ids | ranges::to<std::vector>,
        dof_table_boundary, variable_id, component_id, t, x, bc_values);
}
}  // namespace ProcessLib
