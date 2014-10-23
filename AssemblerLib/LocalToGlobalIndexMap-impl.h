/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

namespace AssemblerLib
{

template <AssemblerLib::ComponentOrder Order>
void
LocalToGlobalIndexMap::constructGlobalIndicesForMeshSubsets()
{
    // For all MeshSubsets and each of their MeshSubset's and each element
    // of that MeshSubset save a line of global indices.
    for (MeshLib::MeshSubsets const* const mss : _mesh_subsets)
    {
        for (MeshLib::MeshSubset const* const ms : *mss)
        {
            std::size_t const mesh_id = ms->getMeshID();

            // For each element find the global indices for node/element
            // components.
            for (auto e = ms->elementsBegin();
                    e != ms->elementsEnd(); ++e)
            {
                appendGlobalIndices<Order>(mesh_id, **e);
            }
        }
    }
}

template <AssemblerLib::ComponentOrder Order>
void
LocalToGlobalIndexMap::appendGlobalIndices(std::size_t const mesh_id, MeshLib::Element const& e)
{
    std::vector<MeshLib::Location> vec_items;
    std::size_t const nnodes = e.getNNodes();
    vec_items.reserve(nnodes);

    for (std::size_t n = 0; n < nnodes; n++)
    {
        vec_items.emplace_back(
            mesh_id,
            MeshLib::MeshItemType::Node,
            e.getNode(n)->getID());
    }

    // Save a line of indices for the current element.
    _rows.push_back(_mesh_component_map.getGlobalIndices<Order>(vec_items));
}

}
