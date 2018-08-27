/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BoundaryCondition.h"

#include "MeshLib/PropertyVector.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
class NonuniformDirichletBoundaryCondition final : public BoundaryCondition
{
public:
    NonuniformDirichletBoundaryCondition(
        // int const bulk_mesh_dimension,
        MeshLib::Mesh const& boundary_mesh,
        MeshLib::PropertyVector<double> const& values,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id,
        int const component_id)
        : _values(values),
          _boundary_mesh(boundary_mesh),
          _variable_id(variable_id),
          _component_id(component_id)
    {
        if (_variable_id >=
                static_cast<int>(dof_table_bulk.getNumberOfVariables()) ||
            _component_id >=
                dof_table_bulk.getNumberOfVariableComponents(_variable_id))
        {
            OGS_FATAL(
                "Variable id or component id too high. Actual values: (%d, "
                "%d), "
                "maximum values: (%d, %d).",
                _variable_id, _component_id,
                dof_table_bulk.getNumberOfVariables(),
                dof_table_bulk.getNumberOfVariableComponents(_variable_id));
        }

        if (!_boundary_mesh.getProperties().existsPropertyVector<std::size_t>(
                "bulk_node_ids"))
        {
            OGS_FATAL(
                "The required bulk node ids map does not exist in the boundary "
                "mesh '%s'.",
                _boundary_mesh.getName().c_str());
        }

        std::vector<MeshLib::Node*> const& bc_nodes = _boundary_mesh.getNodes();
        DBUG(
            "Found %d nodes for Natural BCs for the variable %d and component "
            "%d",
            bc_nodes.size(), variable_id, component_id);

        MeshLib::MeshSubset boundary_mesh_subset(_boundary_mesh, bc_nodes);

        // Create local DOF table from the BC mesh subset for the given variable
        // and component id.
        _dof_table_boundary.reset(dof_table_bulk.deriveBoundaryConstrainedMap(
            variable_id, {component_id}, std::move(boundary_mesh_subset)));
    }

    void getEssentialBCValues(
        const double /*t*/, GlobalVector const& /*x*/,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override
    {
        bc_values.ids.clear();
        bc_values.values.clear();

        // Convert mesh node ids to global index for the given component.
        bc_values.ids.reserve(_values.size());
        bc_values.values.reserve(_values.size());

        // Map boundary dof indices to bulk dof indices and the corresponding
        // values.
        for (auto const* const node : _boundary_mesh.getNodes())
        {
            auto const node_id = node->getID();
            auto const global_index = _dof_table_boundary->getGlobalIndex(
                {_boundary_mesh.getID(), MeshLib::MeshItemType::Node, node_id},
                _variable_id, _component_id);
            if (global_index == NumLib::MeshComponentMap::nop)
                continue;
            // For the DDC approach (e.g. with PETSc option), the negative index
            // of global_index means that the entry by that index is a ghost
            // one, which should be dropped. Especially for PETSc routines
            // MatZeroRows and MatZeroRowsColumns, which are called to apply the
            // Dirichlet BC, the negative index is not accepted like other
            // matrix or vector PETSc routines. Therefore, the following
            // if-condition is applied.
            if (global_index >= 0)
            {
                bc_values.ids.emplace_back(global_index);
                bc_values.values.push_back(_values[node_id]);
            }
        }
    }

private:
    MeshLib::PropertyVector<double> const& _values;
    MeshLib::Mesh const& _boundary_mesh;
    std::unique_ptr<NumLib::LocalToGlobalIndexMap const> _dof_table_boundary;
    int const _variable_id;
    int const _component_id;
};

std::unique_ptr<NonuniformDirichletBoundaryCondition>
createNonuniformDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& boundary_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, const MeshLib::Mesh& bulk_mesh);

}  // ProcessLib
