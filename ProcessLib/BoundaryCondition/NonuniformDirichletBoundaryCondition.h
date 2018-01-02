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
        std::unique_ptr<MeshLib::Mesh>
            boundary_mesh,
        MeshLib::PropertyVector<double> const& values,
        std::size_t const bulk_mesh_id,
        MeshLib::PropertyVector<std::size_t> const& mapping_to_bulk_nodes,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id_bulk,
        int const component_id_bulk)
        : _bulk_mesh_id(bulk_mesh_id),
          _values(values),
          _boundary_mesh(std::move(boundary_mesh)),
          _mapping_to_bulk_nodes(mapping_to_bulk_nodes),
          _dof_table_bulk(dof_table_bulk),
          _variable_id_bulk(variable_id_bulk),
          _component_id_bulk(component_id_bulk)
    {
        if (_variable_id_bulk >=
                static_cast<int>(_dof_table_bulk.getNumberOfVariables()) ||
            _component_id_bulk >= _dof_table_bulk.getNumberOfVariableComponents(
                                      _variable_id_bulk))
        {
            OGS_FATAL(
                "Variable id or component id too high. Actual values: (%d, "
                "%d), "
                "maximum values: (%d, %d).",
                _variable_id_bulk, _component_id_bulk,
                _dof_table_bulk.getNumberOfVariables(),
                _dof_table_bulk.getNumberOfVariableComponents(
                    _variable_id_bulk));
        }
    }

    void getEssentialBCValues(
        const double /*t*/,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override
    {
        SpatialPosition pos;

        bc_values.ids.clear();
        bc_values.values.clear();

        // Convert mesh node ids to global index for the given component.
        bc_values.ids.reserve(_values.size());
        bc_values.values.reserve(_values.size());

        // Map boundary dof indices to bulk dof indices and the corresponding
        // values.
        for (std::size_t i = 0; i < _boundary_mesh->getNumberOfNodes(); ++i)
        {
            auto const bulk_node_id = _mapping_to_bulk_nodes.getComponent(i, 0);

            MeshLib::Location const l{
                _bulk_mesh_id, MeshLib::MeshItemType::Node, bulk_node_id};

            auto const global_index = _dof_table_bulk.getGlobalIndex(
                l, _variable_id_bulk, _component_id_bulk);
            assert(global_index != NumLib::MeshComponentMap::nop);

            bc_values.ids.push_back(global_index);
            bc_values.values.push_back(_values.getComponent(i, 0));
        }
    }

private:
    std::size_t _bulk_mesh_id;
    MeshLib::PropertyVector<double> const& _values;
    std::unique_ptr<MeshLib::Mesh> _boundary_mesh;
    MeshLib::PropertyVector<std::size_t> const& _mapping_to_bulk_nodes;
    NumLib::LocalToGlobalIndexMap const& _dof_table_bulk;
    int const _variable_id_bulk;
    int const _component_id_bulk;
};

std::unique_ptr<NonuniformDirichletBoundaryCondition>
createNonuniformDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, const MeshLib::Mesh& bulk_mesh);

}  // ProcessLib
