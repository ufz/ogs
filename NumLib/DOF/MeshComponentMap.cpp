/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MeshComponentMap.h"

#include "BaseLib/Error.h"
#include "MeshLib/MeshSubset.h"
#include "MeshLib/Node.h"

#ifdef USE_PETSC
#include "MeshLib/NodePartitionedMesh.h"
#endif

namespace NumLib
{
using namespace detail;

#ifdef USE_PETSC
MeshComponentMap::MeshComponentMap(
    std::vector<MeshLib::MeshSubset> const& components, ComponentOrder order)
{
    // Use PETSc with single thread
    const MeshLib::NodePartitionedMesh& partitioned_mesh =
        static_cast<const MeshLib::NodePartitionedMesh&>(
            components[0].getMesh());
    if (partitioned_mesh.isForSingleThread())
    {
        createSerialMeshComponentMap(components, order);
        return;
    }
    else
    {
        createParallelMeshComponentMap(components, order);
    }
}
#else
MeshComponentMap::MeshComponentMap(
    std::vector<MeshLib::MeshSubset> const& components, ComponentOrder order)
{
    createSerialMeshComponentMap(components, order);
}
#endif  // end of USE_PETSC

MeshComponentMap MeshComponentMap::getSubset(
    std::vector<MeshLib::MeshSubset> const& bulk_mesh_subsets,
    MeshLib::MeshSubset const& new_mesh_subset,
    std::vector<int> const& new_global_component_ids) const
{
    {  // Testing first an assumption met later in the code that the meshes for
       // all the bulk_mesh_subsets are equal.
        auto const first_mismatch =
            std::adjacent_find(begin(bulk_mesh_subsets), end(bulk_mesh_subsets),
                               [](auto const& a, auto const& b)
                               { return a.getMeshID() != b.getMeshID(); });
        if (first_mismatch != end(bulk_mesh_subsets))
        {
            OGS_FATAL(
                "Assumption in the MeshComponentMap violated. Expecting all of "
                "mesh ids to be the same, but it is not true for "
                "the mesh '{:s}' with id {:d}.",
                first_mismatch->getMesh().getName(),
                first_mismatch->getMeshID());
        }
    }

    // Mapping of the nodes in the new_mesh_subset to the bulk mesh nodes
    auto bulk_node_ids = [](auto const& mesh)
    {
        auto const* bulk_node_ids_ptr = MeshLib::bulkNodeIDs(mesh);
        if (bulk_node_ids_ptr == nullptr)
        {
            OGS_FATAL(
                "Bulk node ids map expected in the construction of the mesh "
                "subset.");
        }
        return *bulk_node_ids_ptr;
    };
    auto const& bulk_node_ids_map = bulk_node_ids(new_mesh_subset.getMesh());

    // New dictionary for the subset.
    ComponentGlobalIndexDict subset_dict;

    std::size_t const new_mesh_id = new_mesh_subset.getMeshID();
    // Lookup the locations in the current mesh component map and
    // insert the full lines into the new subset dictionary.
    for (auto* const node : new_mesh_subset.getNodes())
    {
        auto const node_id = node->getID();
        bool const is_base_node = MeshLib::isBaseNode(
            *node,
            new_mesh_subset.getMesh().getElementsConnectedToNode(node_id));

        MeshLib::Location const new_location{
            new_mesh_id, MeshLib::MeshItemType::Node, node_id};

        // Assuming the meshes for all the bulk_mesh_subsets are equal.
        MeshLib::Location const bulk_location{
            bulk_mesh_subsets.front().getMeshID(), MeshLib::MeshItemType::Node,
            bulk_node_ids_map[node_id]};

        for (auto component_id : new_global_component_ids)
        {
            auto const global_index =
                getGlobalIndex(bulk_location, component_id);
            if (global_index == nop)
            {
                if (is_base_node)
                {
                    OGS_FATAL(
                        "Could not find a global index for global component "
                        "{:d} for the mesh '{:s}', node {:d}, in the "
                        "corresponding bulk mesh '{:s}' and node {:d}. This "
                        "happens because the boundary mesh is larger then the "
                        "definition region of the bulk component, usually "
                        "because the geometry for the boundary condition is "
                        "too large.",
                        component_id,
                        new_mesh_subset.getMesh().getName(),
                        node_id,
                        bulk_mesh_subsets.front().getMesh().getName(),
                        bulk_node_ids_map[node_id]);
                }
                continue;
            }
            subset_dict.insert({new_location, component_id, global_index});
        }
    }

    return MeshComponentMap(subset_dict);
}

void MeshComponentMap::renumberByLocation(GlobalIndexType offset)
{
    GlobalIndexType global_index = offset;

    auto& m = _dict.get<ByLocation>();  // view as sorted by mesh item
    for (auto itr_mesh_item = m.begin(); itr_mesh_item != m.end();
         ++itr_mesh_item)
    {
        Line pos = *itr_mesh_item;
        pos.global_index = global_index++;
        m.replace(itr_mesh_item, pos);
    }
}

std::vector<int> MeshComponentMap::getComponentIDs(const Location& l) const
{
    auto const& m = _dict.get<ByLocation>();
    auto const p = m.equal_range(Line(l));
    std::vector<int> vec_compID;
    for (auto itr = p.first; itr != p.second; ++itr)
    {
        vec_compID.push_back(itr->comp_id);
    }
    return vec_compID;
}

GlobalIndexType MeshComponentMap::getGlobalIndex(Location const& l,
                                                 int const comp_id) const
{
    auto const& m = _dict.get<ByLocationAndComponent>();
    auto const itr = m.find(Line(l, comp_id));
    return itr != m.end() ? itr->global_index : nop;
}

std::vector<GlobalIndexType> MeshComponentMap::getGlobalIndices(
    const Location& l) const
{
    auto const& m = _dict.get<ByLocation>();
    auto const p = m.equal_range(Line(l));
    std::vector<GlobalIndexType> global_indices;
    for (auto itr = p.first; itr != p.second; ++itr)
    {
        global_indices.push_back(itr->global_index);
    }
    return global_indices;
}

std::vector<GlobalIndexType> MeshComponentMap::getGlobalIndicesByLocation(
    std::vector<Location> const& ls) const
{
    // Create vector of global indices sorted by location containing all
    // locations given in ls parameter.

    std::vector<GlobalIndexType> global_indices;
    global_indices.reserve(ls.size());

    auto const& m = _dict.get<ByLocation>();
    for (const auto& l : ls)
    {
        auto const p = m.equal_range(Line(l));
        for (auto itr = p.first; itr != p.second; ++itr)
        {
            global_indices.push_back(itr->global_index);
        }
    }

    return global_indices;
}

std::vector<GlobalIndexType> MeshComponentMap::getGlobalIndicesByComponent(
    std::vector<Location> const& ls) const
{
    // vector of (Component, global Index) pairs.
    using CIPair = std::pair<int, GlobalIndexType>;
    std::vector<CIPair> pairs;
    pairs.reserve(ls.size());

    // Create a sub dictionary containing all lines with location from ls.
    auto const& m = _dict.get<ByLocation>();
    for (const auto& l : ls)
    {
        auto const p = m.equal_range(Line(l));
        for (auto itr = p.first; itr != p.second; ++itr)
        {
            pairs.emplace_back(itr->comp_id, itr->global_index);
        }
    }

    auto CIPairLess = [](CIPair const& a, CIPair const& b)
    { return a.first < b.first; };

    // Create vector of global indices from sub dictionary sorting by component.
    if (!std::is_sorted(pairs.begin(), pairs.end(), CIPairLess))
    {
        std::stable_sort(pairs.begin(), pairs.end(), CIPairLess);
    }

    std::vector<GlobalIndexType> global_indices;
    global_indices.reserve(pairs.size());
    transform(cbegin(pairs), cend(pairs), back_inserter(global_indices),
              [&](const auto& pair) { return pair.second; });

    return global_indices;
}

GlobalIndexType MeshComponentMap::getLocalIndex(
    Location const& l,
    int const comp_id,
    std::size_t const range_begin,
    std::size_t const range_end) const
{
    GlobalIndexType const global_index = getGlobalIndex(l, comp_id);
    // request for index of linear quantities at higher order nodes
    // results in returning nop
    // That index shall not be modified like a usual global index.
    if (global_index == nop)
    {
        return nop;
    }
#ifndef USE_PETSC
    (void)range_begin;
    (void)range_end;
    return global_index;
#else
    if (global_index >= 0)  // non-ghost location.
        return global_index - range_begin;

    //
    // For a ghost location look up the global index in ghost indices.
    //

    // A special case for a ghost location with global index equal to the size
    // of the local vector:
    GlobalIndexType const real_global_index =
        (-global_index == static_cast<GlobalIndexType>(_num_global_dof))
            ? 0
            : -global_index;

    // TODO Find in ghost indices is O(n^2/2) for n being the length of
    // _ghosts_indices. Providing an inverted table would be faster.
    auto const ghost_index_it = std::find(
        _ghosts_indices.begin(), _ghosts_indices.end(), real_global_index);
    if (ghost_index_it == _ghosts_indices.end())
    {
        OGS_FATAL("index {:d} not found in ghost_indices", real_global_index);
    }

    // Using std::distance on a std::vector is O(1). As long as _ghost_indices
    // remains of std::vector type, this shall be fast.
    return range_end - range_begin +
           std::distance(_ghosts_indices.begin(), ghost_index_it);

#endif
}

void MeshComponentMap::createSerialMeshComponentMap(
    std::vector<MeshLib::MeshSubset> const& components, ComponentOrder order)
{
    // construct dict (and here we number global_index by component type)
    GlobalIndexType global_index = 0;
    int comp_id = 0;
    for (auto const& c : components)
    {
        std::size_t const mesh_id = c.getMeshID();
        auto const& mesh_subset_nodes = c.getNodes();
        // mesh items are ordered first by node, cell, ....
        for (std::size_t j = 0; j < mesh_subset_nodes.size(); j++)
        {
            auto const node_id = mesh_subset_nodes[j]->getID();
            _dict.insert(
                Line(Location(mesh_id, MeshLib::MeshItemType::Node, node_id),
                     comp_id, global_index++));
        }
        comp_id++;
    }
    _num_local_dof = _dict.size();

    if (order == ComponentOrder::BY_LOCATION)
    {
        renumberByLocation();
    }
}

#ifdef USE_PETSC

GlobalIndexType getGlobalIndexWithTaylorHoodElement(
    MeshLib::NodePartitionedMesh const& partitioned_mesh,
    std::size_t const global_node_id,
    int const number_of_components_at_base_node,
    int const number_of_components_at_high_order_node,
    int const component_id_at_base_node,
    int const component_id_at_high_order_node,
    bool const is_base_node)
{
    int const partition_id = partitioned_mesh.getPartitionID(global_node_id);

    auto const n_total_active_base_nodes_before_this_rank =
        partitioned_mesh.getNumberOfRegularBaseNodesAtRank(partition_id);

    auto const n_total_active_high_order_nodes_before_this_rank =
        partitioned_mesh.getNumberOfRegularHighOrderNodesAtRank(partition_id);

    auto const node_id_offset =
        n_total_active_base_nodes_before_this_rank +
        n_total_active_high_order_nodes_before_this_rank;

    auto const index_offset = n_total_active_base_nodes_before_this_rank *
                                  number_of_components_at_base_node +
                              n_total_active_high_order_nodes_before_this_rank *
                                  number_of_components_at_high_order_node;

    if (is_base_node)
    {
        return static_cast<GlobalIndexType>(
                   index_offset + (global_node_id - node_id_offset) *
                                      number_of_components_at_base_node) +
               component_id_at_base_node;
    }

    int const n_active_base_nodes_of_this_partition =
        partitioned_mesh.getNumberOfRegularBaseNodesAtRank(partition_id + 1) -
        n_total_active_base_nodes_before_this_rank;

    /*
        The global indices of components are numbered as what depicted below by
        assuming that the base node has three components and the high order node
        has two components:

        Partition     |       0       |      1       |   ...
                      --------------------------------
        Regular nodes | Base | higher | Base | higher|   ...
                      --------------------------------
                  c0  x      x         x      x          ...
                  c1  x      x         x      x          ...
                  c2  x                x                 ...
    */
    return static_cast<GlobalIndexType>(
               index_offset +
               n_active_base_nodes_of_this_partition *
                   number_of_components_at_base_node +
               (global_node_id - node_id_offset -
                n_active_base_nodes_of_this_partition) *
                   number_of_components_at_high_order_node) +
           component_id_at_high_order_node;
}

void MeshComponentMap::createParallelMeshComponentMap(
    std::vector<MeshLib::MeshSubset> const& components, ComponentOrder order)
{
    if (order != ComponentOrder::BY_LOCATION)
    {
        // Not allowed in parallel case since this is not suitable to
        // arrange non ghost entries of a partition within
        // a rank in the parallel computing.
        OGS_FATAL(
            "Global index in the system of equations can only be numbered by "
            "the order type of ComponentOrder::BY_LOCATION");
    }

    const MeshLib::NodePartitionedMesh& partitioned_mesh =
        static_cast<const MeshLib::NodePartitionedMesh&>(
            components[0].getMesh());

    //
    // get the number of unknowns and the number of components at extra node
    //
    GlobalIndexType num_unknowns = 0;
    int components_at_high_order_nodes = 0;
    for (auto const& c : components)
    {
        if (partitioned_mesh.getNumberOfNodes() == c.getNodes().size())
        {
            components_at_high_order_nodes++;
            num_unknowns += partitioned_mesh.getNumberOfGlobalNodes();
        }
        else
        {
            num_unknowns += partitioned_mesh.getNumberOfGlobalBaseNodes();
        }
    }

    // construct dict (and here we number global_index by component type)
    //
    int const n_components = components.size();

    int comp_id = 0;
    int comp_id_at_high_order_node = 0;
    _num_global_dof = 0;
    _num_local_dof = 0;
    for (auto const& c : components)
    {
        assert(dynamic_cast<MeshLib::NodePartitionedMesh const*>(
                   &c.getMesh()) != nullptr);
        std::size_t const mesh_id = c.getMeshID();
        const MeshLib::NodePartitionedMesh& partitioned_mesh =
            static_cast<const MeshLib::NodePartitionedMesh&>(c.getMesh());
        const auto& sub_mesh_nodes = c.getNodes();

        // mesh items are ordered first by node, cell, ....
        for (std::size_t j = 0; j < sub_mesh_nodes.size(); j++)
        {
            const auto node_id = sub_mesh_nodes[j]->getID();
            const bool is_base_node =
                MeshLib::isBaseNode(*sub_mesh_nodes[j],
                                    partitioned_mesh.getElementsConnectedToNode(
                                        *sub_mesh_nodes[j]));

            const auto global_node_id =
                partitioned_mesh.getGlobalNodeID(node_id);
            GlobalIndexType global_index =
                (!c.useTaylorHoodElements())
                    ? static_cast<GlobalIndexType>(
                          n_components * global_node_id + comp_id)
                    : getGlobalIndexWithTaylorHoodElement(
                          partitioned_mesh, global_node_id, n_components,
                          components_at_high_order_nodes, comp_id,
                          comp_id_at_high_order_node, is_base_node);
            const bool is_ghost = partitioned_mesh.isGhostNode(node_id);
            if (is_ghost)
            {
                _ghosts_indices.push_back(global_index);
                global_index = -global_index;
                // If the ghost entry has an index of 0,
                // its index is set to the negative value of unknowns.
                if (global_index == 0)
                {
                    global_index = -num_unknowns;
                }
            }
            else
            {
                _num_local_dof++;
            }

            _dict.insert(
                Line(Location(mesh_id, MeshLib::MeshItemType::Node, node_id),
                     comp_id, global_index));
        }

        bool const use_whole_nodes =
            (partitioned_mesh.getNumberOfNodes() == c.getNodes().size());
        if (use_whole_nodes)
        {
            _num_global_dof += partitioned_mesh.getNumberOfGlobalNodes();
            comp_id_at_high_order_node++;
        }
        else
        {
            _num_global_dof += partitioned_mesh.getNumberOfGlobalBaseNodes();
        }

        comp_id++;
    }
}
#endif

}  // namespace NumLib
