/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "BoundaryCondition.h"

namespace ProcessLib
{
// TODO docu
/// The DirichletBoundaryCondition class describes a constant in space
/// and time Dirichlet boundary condition.
/// The expected parameter in the passed configuration is "value" which, when
/// not present defaults to zero.
class DirichletBoundaryCondition final : public BoundaryCondition
{
public:
    DirichletBoundaryCondition(
        Parameter<double> const& parameter, MeshLib::Mesh const& bc_mesh,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id)
        : _parameter(parameter),
          _bc_mesh(bc_mesh),
          _variable_id(variable_id),
          _component_id(component_id)
    {
        if (variable_id >=
                static_cast<int>(dof_table_bulk.getNumberOfVariables()) ||
            component_id >=
                dof_table_bulk.getNumberOfVariableComponents(variable_id))
        {
            OGS_FATAL(
                "Variable id or component id too high. Actual values: (%d, "
                "%d), maximum values: (%d, %d).",
                variable_id, component_id,
                dof_table_bulk.getNumberOfVariables(),
                dof_table_bulk.getNumberOfVariableComponents(variable_id));
        }

        if (!_bc_mesh.getProperties().existsPropertyVector<std::size_t>(
                "bulk_node_ids"))
        {
            OGS_FATAL(
                "The required bulk node ids map does not exist in the boundary "
                "mesh '%s'.", _bc_mesh.getName().c_str());
        }

        std::vector<MeshLib::Node*> const& bc_nodes = _bc_mesh.getNodes();
        DBUG(
            "Found %d nodes for Dirichlet BCs for the variable %d and "
            "component "
            "%d",
            bc_nodes.size(), variable_id, component_id);

        MeshLib::MeshSubset bc_mesh_subset(_bc_mesh, bc_nodes);

        // Create local DOF table from the bc mesh subset for the given variable
        // and component id.
        _dof_table_boundary.reset(dof_table_bulk.deriveBoundaryConstrainedMap(
            variable_id, {component_id}, std::move(bc_mesh_subset)));
    }

    void getEssentialBCValues(
        const double t, GlobalVector const& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

private:
    Parameter<double> const& _parameter;

    MeshLib::Mesh const& _bc_mesh;
    std::unique_ptr<NumLib::LocalToGlobalIndexMap const> _dof_table_boundary;
    int const _variable_id;
    int const _component_id;
};

std::unique_ptr<DirichletBoundaryCondition> createDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    int const component_id,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters);

}  // namespace ProcessLib
