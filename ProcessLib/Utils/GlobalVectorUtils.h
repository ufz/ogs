/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshLib/Properties.h"

#include "NumLib/DOF/LocalToGlobalIndexMap.h"

namespace ProcessLib
{
/// Copies part of a global vector for the given variable into output vector
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
    std::fill(begin(output_vector), end(output_vector),
              std::numeric_limits<double>::quiet_NaN());

    int const n_components =
        local_to_global_index_map.getNumberOfVariableComponents(variable_id);
    for (int component = 0; component < n_components; ++component)
    {
        auto const& mesh_subsets =
            local_to_global_index_map.getMeshSubsets(variable_id, component);
        assert(mesh_subsets.size() ==
               1);  // Multiple meshes are not supported by the output_vector.
        for (auto const& ms : mesh_subsets)
        {
            auto const mesh_id = ms->getMeshID();
            for (auto const& node : ms->getNodes())
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
}
}  // namespace ProcessLib
