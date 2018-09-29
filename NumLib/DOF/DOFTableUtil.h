/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MathLib/LinAlg/LinAlg.h"
#include "MathLib/LinAlg/LinAlgEnums.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"

namespace NumLib
{
//! Returns the value for the given \c node_id and \c global_component_id from
//! the vector \c x in case the node is not a ghost node. Else 0.0 is returned.
double getNonGhostNodalValue(GlobalVector const& x, MeshLib::Mesh const& mesh,
                             NumLib::LocalToGlobalIndexMap const& dof_table,
                             std::size_t const node_id,
                             std::size_t const global_component_id);

//! Returns the value for the given \c node_id and \c global_component_id from
//! the vector \c x.
double getNodalValue(GlobalVector const& x, MeshLib::Mesh const& mesh,
                     NumLib::LocalToGlobalIndexMap const& dof_table,
                     std::size_t const node_id,
                     std::size_t const global_component_id);

//! Returns nodal indices for the item identified by \c mesh_item_id from the
//! given \c dof_table.
std::vector<GlobalIndexType> getIndices(
    std::size_t const mesh_item_id,
    NumLib::LocalToGlobalIndexMap const& dof_table);

//! Returns row/column indices for the item identified by \c id from the
//! given \c dof_table.
LocalToGlobalIndexMap::RowColumnIndices getRowColumnIndices(
    std::size_t const id,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<GlobalIndexType>& indices);

//! Computes the specified norm of the given global component of the given vector x.
//! \remark
//! \c x is typically the solution vector of a monolithically coupled process
//! with several primary variables.
double norm(GlobalVector const& x, unsigned const global_component,
            MathLib::VecNormType norm_type,
            LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh);

/// Copies part of a global vector for the given variable into output_vector
/// while applying a function to each value.
///
/// \attention The output_vector is accessed through node id and component,
/// therefore multiple meshes are not supported.
template <typename Functor>
void transformVariableFromGlobalVector(
    GlobalVector const& input_vector, int const variable_id,
    NumLib::LocalToGlobalIndexMap const& local_to_global_index_map,
    MeshLib::PropertyVector<double>& output_vector, Functor mapFunction)
{
    MathLib::LinAlg::setLocalAccessibleVector(input_vector);

    std::fill(begin(output_vector), end(output_vector),
              std::numeric_limits<double>::quiet_NaN());

    int const n_components =
        local_to_global_index_map.getNumberOfVariableComponents(variable_id);
    for (int component = 0; component < n_components; ++component)
    {
        auto const& mesh_subset =
            local_to_global_index_map.getMeshSubset(variable_id, component);
        auto const mesh_id = mesh_subset.getMeshID();
        for (auto const& node : mesh_subset.getNodes())
        {
            auto const node_id = node->getID();
            MeshLib::Location const l(mesh_id, MeshLib::MeshItemType::Node,
                                      node_id);
            output_vector.getComponent(node_id, component) = mapFunction(
                input_vector[local_to_global_index_map.getGlobalIndex(
                    l, variable_id, component)]);
        }
    }
}
}  // namespace NumLib
