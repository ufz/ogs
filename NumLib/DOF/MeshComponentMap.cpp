/**
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshComponentMap.h"

#include "MeshLib/MeshSubsets.h"

#ifdef USE_PETSC
#include "MeshLib/NodePartitionedMesh.h"
#endif

namespace NumLib
{

using namespace detail;

GlobalIndexType const MeshComponentMap::nop =
    std::numeric_limits<GlobalIndexType>::max();

#ifdef USE_PETSC
MeshComponentMap::MeshComponentMap(
    const std::vector<std::unique_ptr<MeshLib::MeshSubsets>>& components,
    ComponentOrder order)
{
    // get number of unknows
    GlobalIndexType num_unknowns = 0;
    for (auto const& c : components)
    {
        assert(c != nullptr);
        for (unsigned mesh_subset_index = 0; mesh_subset_index < c->size();
             mesh_subset_index++)
        {
            MeshLib::MeshSubset const& mesh_subset =
                c->getMeshSubset(mesh_subset_index);
            // PETSc always works with MeshLib::NodePartitionedMesh.
            const MeshLib::NodePartitionedMesh& mesh =
                static_cast<const MeshLib::NodePartitionedMesh&>(
                    mesh_subset.getMesh());
            num_unknowns += mesh.getNumberOfGlobalNodes();
        }
    }

    // construct dict (and here we number global_index by component type)
    std::size_t cell_index = 0;
    std::size_t comp_id = 0;
    _num_global_dof = 0;
    _num_local_dof = 0;
    for (auto const& c : components)
    {
        assert(c != nullptr);
        for (unsigned mesh_subset_index = 0; mesh_subset_index < c->size();
             mesh_subset_index++)
        {
            MeshLib::MeshSubset const& mesh_subset =
                c->getMeshSubset(mesh_subset_index);
            assert(dynamic_cast<MeshLib::NodePartitionedMesh const*>(
                       &mesh_subset.getMesh()) != nullptr);
            std::size_t const mesh_id = mesh_subset.getMeshID();
            const MeshLib::NodePartitionedMesh& mesh =
                static_cast<const MeshLib::NodePartitionedMesh&>(
                    mesh_subset.getMesh());

            // mesh items are ordered first by node, cell, ....
            for (std::size_t j = 0; j < mesh_subset.getNumberOfNodes(); j++)
            {
                GlobalIndexType global_id = 0;
                if (order == ComponentOrder::BY_LOCATION)
                {
                    global_id = static_cast<GlobalIndexType>(
                        components.size() * mesh.getGlobalNodeID(j) + comp_id);
                }
                else
                {
                    // _num_global_dof is used as the global index offset
                    global_id = static_cast<GlobalIndexType>(
                        _num_global_dof + mesh.getGlobalNodeID(j));
                }
                const bool is_ghost =
                    mesh.isGhostNode(mesh.getNode(j)->getID());
                if (is_ghost)
                {
                    _ghosts_indices.push_back(global_id);
                    global_id = -global_id;
                    // If the ghost entry has an index of 0,
                    // its index is set to the negative value of unknowns.
                    if (global_id == 0) global_id = -num_unknowns;
                }
                else
                    _num_local_dof++;

                _dict.insert(
                    Line(Location(mesh_id, MeshLib::MeshItemType::Node, j),
                         comp_id, global_id));
            }

            // Note: If the cells are really used (e.g. for the mixed FEM),
            // the following global cell index must be reconsidered
            // according to the employed cell indexing method.
            for (std::size_t j = 0; j < mesh_subset.getNumberOfElements(); j++)
                _dict.insert(
                    Line(Location(mesh_id, MeshLib::MeshItemType::Cell, j),
                         comp_id, cell_index++));

            _num_global_dof += mesh.getNumberOfGlobalNodes();
        }
        comp_id++;
    }
}
#else
MeshComponentMap::MeshComponentMap(
    const std::vector<std::unique_ptr<MeshLib::MeshSubsets>>& components,
    ComponentOrder order)
{
    // construct dict (and here we number global_index by component type)
    GlobalIndexType global_index = 0;
    std::size_t comp_id = 0;
    for (auto const& c : components)
    {
        assert (c != nullptr);
        for (std::size_t mesh_subset_index = 0; mesh_subset_index < c->size(); mesh_subset_index++)
        {
            MeshLib::MeshSubset const& mesh_subset = c->getMeshSubset(mesh_subset_index);
            std::size_t const mesh_id = mesh_subset.getMeshID();
            // mesh items are ordered first by node, cell, ....
            for (std::size_t j=0; j<mesh_subset.getNumberOfNodes(); j++)
                _dict.insert(Line(Location(mesh_id, MeshLib::MeshItemType::Node, mesh_subset.getNodeID(j)), comp_id, global_index++));
            for (std::size_t j=0; j<mesh_subset.getNumberOfElements(); j++)
                _dict.insert(Line(Location(mesh_id, MeshLib::MeshItemType::Cell, mesh_subset.getElementID(j)), comp_id, global_index++));
        }
        comp_id++;
    }
    _num_local_dof = _dict.size();

    if (order == ComponentOrder::BY_LOCATION)
        renumberByLocation();
}
#endif // end of USE_PETSC

MeshComponentMap MeshComponentMap::getSubset(
    std::vector<int> const& component_ids,
    MeshLib::MeshSubsets const& mesh_subsets) const
{
    // New dictionary for the subset.
    ComponentGlobalIndexDict subset_dict;

    for (auto const& mesh_subset : mesh_subsets)
    {
        std::size_t const mesh_id = mesh_subset->getMeshID();
        // Lookup the locations in the current mesh component map and
        // insert the full lines into the subset dictionary.
        for (std::size_t j = 0; j < mesh_subset->getNumberOfNodes(); j++)
            for (auto component_id : component_ids)
                subset_dict.insert(
                    getLine(Location(mesh_id, MeshLib::MeshItemType::Node,
                                     mesh_subset->getNodeID(j)),
                            component_id));
        for (std::size_t j = 0; j < mesh_subset->getNumberOfElements(); j++)
            for (auto component_id : component_ids)
                subset_dict.insert(
                    getLine(Location(mesh_id, MeshLib::MeshItemType::Cell,
                                     mesh_subset->getElementID(j)),
                            component_id));
    }

    return MeshComponentMap(subset_dict);
}

void MeshComponentMap::renumberByLocation(GlobalIndexType offset)
{
    GlobalIndexType global_index = offset;

    auto &m = _dict.get<ByLocation>(); // view as sorted by mesh item
    for (auto itr_mesh_item=m.begin(); itr_mesh_item!=m.end(); ++itr_mesh_item)
    {
        Line pos = *itr_mesh_item;
        pos.global_index = global_index++;
        m.replace(itr_mesh_item, pos);
    }
}

std::vector<std::size_t> MeshComponentMap::getComponentIDs(const Location &l) const
{
    auto const &m = _dict.get<ByLocation>();
    auto const p = m.equal_range(Line(l));
    std::vector<std::size_t> vec_compID;
    for (auto itr=p.first; itr!=p.second; ++itr)
        vec_compID.push_back(itr->comp_id);
    return vec_compID;
}

Line MeshComponentMap::getLine(Location const& l,
    std::size_t const comp_id) const
{
    auto const &m = _dict.get<ByLocationAndComponent>();
    auto const itr = m.find(Line(l, comp_id));
    assert(itr != m.end());     // The line must exist in the current dictionary.
    return *itr;
}

GlobalIndexType MeshComponentMap::getGlobalIndex(Location const& l,
    std::size_t const comp_id) const
{
    auto const &m = _dict.get<ByLocationAndComponent>();
    auto const itr = m.find(Line(l, comp_id));
    return itr!=m.end() ? itr->global_index : nop;
}

std::vector<GlobalIndexType> MeshComponentMap::getGlobalIndices(const Location &l) const
{
    auto const &m = _dict.get<ByLocation>();
    auto const p = m.equal_range(Line(l));
    std::vector<GlobalIndexType> global_indices;
    for (auto itr=p.first; itr!=p.second; ++itr)
        global_indices.push_back(itr->global_index);
    return global_indices;
}

std::vector<GlobalIndexType> MeshComponentMap::getGlobalIndicesByLocation(
    std::vector<Location> const& ls) const
{
    // Create vector of global indices sorted by location containing all
    // locations given in ls parameter.

    std::vector<GlobalIndexType> global_indices;
    global_indices.reserve(ls.size());

    auto const &m = _dict.get<ByLocation>();
    for (const auto& l : ls)
    {
        auto const p = m.equal_range(Line(l));
        for (auto itr = p.first; itr != p.second; ++itr)
            global_indices.push_back(itr->global_index);
    }

    return global_indices;
}

std::vector<GlobalIndexType> MeshComponentMap::getGlobalIndicesByComponent(
    std::vector<Location> const& ls) const
{
    // vector of (Component, global Index) pairs.
    typedef std::pair<std::size_t, GlobalIndexType> CIPair;
    std::vector<CIPair> pairs;
    pairs.reserve(ls.size());

    // Create a sub dictionary containing all lines with location from ls.
    auto const &m = _dict.get<ByLocation>();
    for (const auto& l : ls)
    {
        auto const p = m.equal_range(Line(l));
        for (auto itr = p.first; itr != p.second; ++itr)
            pairs.emplace_back(itr->comp_id, itr->global_index);
    }

    auto CIPairLess = [](CIPair const& a, CIPair const& b)
        {
            return a.first < b.first;
        };

    // Create vector of global indices from sub dictionary sorting by component.
    if (!std::is_sorted(pairs.begin(), pairs.end(), CIPairLess))
        std::stable_sort(pairs.begin(), pairs.end(), CIPairLess);

    std::vector<GlobalIndexType> global_indices;
    global_indices.reserve(pairs.size());
    for (const auto& pair : pairs)
        global_indices.push_back(pair.second);

    return global_indices;
}

GlobalIndexType MeshComponentMap::getLocalIndex(
    Location const& l,
    std::size_t const comp_id,
    std::size_t const range_begin,
    std::size_t const range_end) const
{
    GlobalIndexType const global_index = getGlobalIndex(l, comp_id);
#ifndef USE_PETSC
    (void)range_begin;
    (void)range_end;
    return global_index;
#else
    if (global_index >= 0)    // non-ghost location.
        return global_index - range_begin;

    //
    // For a ghost location look up the global index in ghost indices.
    //

    // A special case for a ghost location with global index equal to the size
    // of the local vector:
    if (-global_index == static_cast<GlobalIndexType>(_num_global_dof))
        return 0;

    // TODO Find in ghost indices is O(n^2/2) for n being the length of
    // _ghosts_indices. Providing an inverted table would be faster.
    auto const ghost_index_it = std::find(_ghosts_indices.begin(),
                                          _ghosts_indices.end(), -global_index);
    if (ghost_index_it == _ghosts_indices.end())
    {
        OGS_FATAL("index %d not found in ghost_indices", -global_index);
    }

    // Using std::distance on a std::vector is O(1). As long as _ghost_indices
    // remains of std::vector type, this shall be fast.
    return range_end - range_begin +
           std::distance(_ghosts_indices.begin(), ghost_index_it);

#endif
}

}   // namespace NumLib
