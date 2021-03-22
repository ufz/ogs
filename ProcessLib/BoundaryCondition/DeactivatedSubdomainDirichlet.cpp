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
#include "ProcessLib/DeactivatedSubdomain.h"

namespace ProcessLib
{
DeactivatedSubdomainDirichlet::DeactivatedSubdomainDirichlet(
    std::vector<std::size_t> const* const active_element_ids,
    BaseLib::TimeInterval const& time_interval,
    ParameterLib::Parameter<double> const& parameter,
    DeactivatedSubdomainMesh const& subdomain,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    int const component_id)
    : _parameter(parameter),
      _subdomain(subdomain),
      _variable_id(variable_id),
      _component_id(component_id),
      _time_interval(time_interval),
      _active_element_ids(active_element_ids)
{
    config(dof_table_bulk);
}

void DeactivatedSubdomainDirichlet::config(
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk)
{
    checkParametersOfDirichletBoundaryCondition(*_subdomain.mesh, dof_table_bulk,
                                                _variable_id, _component_id);

    std::vector<MeshLib::Node*> const& bc_nodes = _subdomain.mesh->getNodes();
    MeshLib::MeshSubset subdomain_mesh_subset(*_subdomain.mesh, bc_nodes);

    // Create local DOF table from the BC mesh subset for the given variable
    // and component id.
    _dof_table_boundary.reset(dof_table_bulk.deriveBoundaryConstrainedMap(
        _variable_id, {_component_id}, std::move(subdomain_mesh_subset)));
}

void DeactivatedSubdomainDirichlet::getEssentialBCValues(
    const double t, GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    auto const& bulk_element_ids =
        *_subdomain.mesh->getProperties()
             .template getPropertyVector<std::size_t>(
                 "bulk_element_ids", MeshLib::MeshItemType::Cell, 1);

    auto is_inactive = [&](MeshLib::Element const* const e) {
        return none_of(
            begin(*_active_element_ids), end(*_active_element_ids),
            [&](auto const id) { return id == bulk_element_ids[e->getID()]; });
    };

    std::vector<MeshLib::Node*> inactive_nodes_in_bc_mesh;
    std::copy_if(begin(_subdomain.inner_nodes), end(_subdomain.inner_nodes),
                 back_inserter(inactive_nodes_in_bc_mesh),
                 [&](MeshLib::Node* const n) {
                     const auto& connected_elements = n->getElements();

                     return std::all_of(begin(connected_elements),
                                        end(connected_elements), is_inactive);
                 });
    if (_time_interval.contains(t))
    {
        getEssentialBCValuesLocal(_parameter, *_subdomain.mesh,
                                  inactive_nodes_in_bc_mesh, *_dof_table_boundary,
                                  _variable_id, _component_id, t, x, bc_values);
        return;
    }

    bc_values.ids.clear();
    bc_values.values.clear();
}
}  // namespace ProcessLib
