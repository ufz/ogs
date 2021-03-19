/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "DeactivatedSubdomainDirichlet.h"

#include "BaseLib/TimeInterval.h"
#include "DirichletBoundaryCondition.h"
#include "DirichletBoundaryConditionAuxiliaryFunctions.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
DeactivatedSubdomainDirichlet::DeactivatedSubdomainDirichlet(
    BaseLib::TimeInterval const& time_interval,
    ParameterLib::Parameter<double> const& parameter,
    MeshLib::Mesh const& bc_mesh,
    std::vector<MeshLib::Node*> const& nodes_in_bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    int const component_id)
    : _parameter(parameter),
      _bc_mesh(bc_mesh),
      _nodes_in_bc_mesh(nodes_in_bc_mesh),
      _variable_id(variable_id),
      _component_id(component_id),
      _time_interval(time_interval)
{
    config(dof_table_bulk);
}

void DeactivatedSubdomainDirichlet::config(
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk)
{
    checkParametersOfDirichletBoundaryCondition(_bc_mesh, dof_table_bulk,
                                                _variable_id, _component_id);

    std::vector<MeshLib::Node*> const& bc_nodes = _bc_mesh.getNodes();
    MeshLib::MeshSubset bc_mesh_subset(_bc_mesh, bc_nodes);

    // Create local DOF table from the BC mesh subset for the given variable
    // and component id.
    _dof_table_boundary.reset(dof_table_bulk.deriveBoundaryConstrainedMap(
        _variable_id, {_component_id}, std::move(bc_mesh_subset)));
}

void DeactivatedSubdomainDirichlet::getEssentialBCValues(
    const double t, GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    if (_time_interval.contains(t))
    {
        getEssentialBCValuesLocal(_parameter, _bc_mesh, _nodes_in_bc_mesh,
                                  *_dof_table_boundary, _variable_id,
                                  _component_id, t, x, bc_values);
        return;
    }

    bc_values.ids.clear();
    bc_values.values.clear();
}
}  // namespace ProcessLib
