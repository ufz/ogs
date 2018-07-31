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
    DirichletBoundaryCondition(Parameter<double> const& parameter,
                               MeshLib::Mesh const& bc_mesh,
                               NumLib::LocalToGlobalIndexMap const& dof_table,
                               std::size_t const bulk_mesh_id,
                               int const variable_id, int const component_id)
        : _parameter(parameter),
          _bc_mesh(bc_mesh),
          _dof_table(dof_table),
          _bulk_mesh_id(bulk_mesh_id),
          _variable_id(variable_id),
          _component_id(component_id)
    {
        if (variable_id >= static_cast<int>(dof_table.getNumberOfVariables()) ||
            component_id >=
                dof_table.getNumberOfVariableComponents(variable_id))
        {
            OGS_FATAL(
                "Variable id or component id too high. Actual values: (%d, "
                "%d), "
                "maximum values: (%d, %d).",
                variable_id, component_id, dof_table.getNumberOfVariables(),
                dof_table.getNumberOfVariableComponents(variable_id));
        }

        if (!_bc_mesh.getProperties().existsPropertyVector<std::size_t>(
                "bulk_node_ids"))
        {
            OGS_FATAL(
                "The required bulk node ids map does not exist in the boundary "
                "mesh '%s'.", _bc_mesh.getName().c_str());
        }
    }

    void getEssentialBCValues(
        const double t, GlobalVector const& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

private:
    Parameter<double> const& _parameter;

    MeshLib::Mesh const& _bc_mesh;
    NumLib::LocalToGlobalIndexMap const& _dof_table;
    std::size_t const _bulk_mesh_id;
    int const _variable_id;
    int const _component_id;
};

std::unique_ptr<DirichletBoundaryCondition> createDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, std::size_t const mesh_id,
    int const variable_id, int const component_id,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters);

}  // namespace ProcessLib
