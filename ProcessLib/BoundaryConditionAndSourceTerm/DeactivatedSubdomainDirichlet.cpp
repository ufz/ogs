/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#include "DeactivatedSubdomainDirichlet.h"

#include "DirichletBoundaryCondition.h"
#include "DirichletBoundaryConditionAuxiliaryFunctions.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/DeactivatedSubdomain.h"

namespace ProcessLib
{
DeactivatedSubdomainDirichlet::DeactivatedSubdomainDirichlet(
    MeshLib::PropertyVector<unsigned char> const& is_active,
    MathLib::PiecewiseLinearInterpolation time_interval,
    ParameterLib::Parameter<double> const& parameter,
    bool const set_outer_nodes_dirichlet_values,
    DeactivatedSubdomainMesh const& subdomain,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    int const component_id)
    : _parameter(parameter),
      _subdomain(subdomain),
      _variable_id(variable_id),
      _component_id(component_id),
      _time_interval(std::move(time_interval)),
      _is_active(is_active),
      _set_outer_nodes_dirichlet_values(set_outer_nodes_dirichlet_values)
{
    config(dof_table_bulk);
}

void DeactivatedSubdomainDirichlet::config(
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk)
{
    checkParametersOfDirichletBoundaryCondition(_subdomain.mesh, dof_table_bulk,
                                                _variable_id, _component_id);

    std::vector<MeshLib::Node*> const& bc_nodes = _subdomain.mesh.getNodes();
    MeshLib::MeshSubset subdomain_mesh_subset(_subdomain.mesh, bc_nodes);

    // Create local DOF table from the BC mesh subset for the given variable
    // and component id.
    _dof_table_boundary = dof_table_bulk.deriveBoundaryConstrainedMap(
        _variable_id, {_component_id}, std::move(subdomain_mesh_subset));
}

void DeactivatedSubdomainDirichlet::getEssentialBCValues(
    const double t, GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    [[maybe_unused]] auto const& bulk_node_ids =
        *MeshLib::bulkNodeIDs(_subdomain.mesh);
    [[maybe_unused]] auto const& bulk_element_ids =
        *MeshLib::bulkElementIDs(_subdomain.mesh);

    auto is_inactive_id = [&](std::size_t const bulk_element_id)
    { return _is_active[bulk_element_id] == 0; };

    auto is_inactive_element = [&](MeshLib::Element const* const e)
    { return is_inactive_id(bulk_element_ids[e->getID()]); };

    std::vector<std::size_t> inactive_nodes_in_bc_mesh;
    std::copy_if(begin(_subdomain.inner_nodes), end(_subdomain.inner_nodes),
                 back_inserter(inactive_nodes_in_bc_mesh),
                 [&](std::size_t const n)
                 {
                     const auto& connected_elements =
                         _subdomain.mesh.getElementsConnectedToNode(n);

                     return std::all_of(begin(connected_elements),
                                        end(connected_elements),
                                        is_inactive_element);
                 });

    if (_set_outer_nodes_dirichlet_values)
    {
        std::copy_if(begin(_subdomain.outer_nodes), end(_subdomain.outer_nodes),
                     back_inserter(inactive_nodes_in_bc_mesh),
                     [&](std::size_t const n)
                     {
                         const auto& connected_elements =
                             _subdomain.mesh.getElementsConnectedToNode(n);

                         return std::all_of(begin(connected_elements),
                                            end(connected_elements),
                                            is_inactive_element);
                     });
    }
    else
    {
        for (std::size_t i = 0; i < _subdomain.outer_nodes.size(); ++i)
        {
            auto const& connected_elements = _subdomain.outer_nodes_elements[i];
            if (std::all_of(begin(connected_elements), end(connected_elements),
                            is_inactive_id))
            {
                inactive_nodes_in_bc_mesh.push_back(_subdomain.outer_nodes[i]);
            }
        }
    }

    auto time_interval_contains = [&](double const t)
    {
        return _time_interval.getSupportMin() <= t &&
               t <= _time_interval.getSupportMax();
    };
    if (time_interval_contains(t))
    {
        getEssentialBCValuesLocal(
            _parameter, _subdomain.mesh, inactive_nodes_in_bc_mesh,
            *_dof_table_boundary, _variable_id, _component_id, t, x, bc_values);
        return;
    }

    bc_values.ids.clear();
    bc_values.values.clear();
}
}  // namespace ProcessLib
