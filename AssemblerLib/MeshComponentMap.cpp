/**
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <iostream>
#include <numeric>

#include <boost/range/algorithm/for_each.hpp>

#include "MeshLib/MeshSubsets.h"

#include "MeshComponentMap.h"

namespace AssemblerLib
{

using namespace detail;

std::size_t const MeshComponentMap::nop = std::numeric_limits<std::size_t>::max();

MeshComponentMap::MeshComponentMap(
    const std::vector<MeshLib::MeshSubsets*> &components, ComponentOrder order)
{
    // construct dict (and here we number global_index by component type)
    std::size_t global_index = 0;
    std::size_t comp_id = 0;
    for (auto c = components.cbegin(); c != components.cend(); ++c)
    {
        for (unsigned mesh_subset_index = 0; mesh_subset_index < (*c)->size(); mesh_subset_index++)
        {
            MeshLib::MeshSubset const& mesh_subset = (*c)->getMeshSubset(mesh_subset_index);
            std::size_t const mesh_id = mesh_subset.getMeshID();
            // mesh items are ordered first by node, cell, ....
            for (std::size_t j=0; j<mesh_subset.getNNodes(); j++)
                _dict.insert(Line(Location(mesh_id, MeshLib::MeshItemType::Node, j), comp_id, global_index++));
            for (std::size_t j=0; j<mesh_subset.getNElements(); j++)
                _dict.insert(Line(Location(mesh_id, MeshLib::MeshItemType::Cell, j), comp_id, global_index++));
        }
        comp_id++;
    }

    if (order == ComponentOrder::BY_LOCATION)
        renumberByLocation();
}

void MeshComponentMap::renumberByLocation(std::size_t offset)
{
    std::size_t global_index = offset;

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

std::size_t MeshComponentMap::getGlobalIndex(Location const& l,
    std::size_t const c) const
{
    auto const &m = _dict.get<ByLocationAndComponent>();
    auto const itr = m.find(Line(l, c));
    return itr!=m.end() ? itr->global_index
                        : std::numeric_limits<std::size_t>::max();
}

std::vector<std::size_t> MeshComponentMap::getGlobalIndices(const Location &l) const
{
    auto const &m = _dict.get<ByLocation>();
    auto const p = m.equal_range(Line(l));
    std::vector<std::size_t> global_indices;
    for (auto itr=p.first; itr!=p.second; ++itr)
        global_indices.push_back(itr->global_index);
    return global_indices;
}

template <>
std::vector<std::size_t>
MeshComponentMap::getGlobalIndices<ComponentOrder::BY_LOCATION>(
    std::vector<Location> const &ls) const
{
    // Create vector of global indices sorted by location containing all
    // locations given in ls parameter.

    std::vector<std::size_t> global_indices;
    global_indices.reserve(ls.size());

    auto const &m = _dict.get<ByLocation>();
    for (auto l = ls.cbegin(); l != ls.cend(); ++l)
    {
        auto const p = m.equal_range(Line(*l));
        for (auto itr = p.first; itr != p.second; ++itr)
            global_indices.push_back(itr->global_index);
    }

    return global_indices;
}

template <>
std::vector<std::size_t>
MeshComponentMap::getGlobalIndices<ComponentOrder::BY_COMPONENT>(
    std::vector<Location> const &ls) const
{
    // vector of (Component, global Index) pairs.
    typedef std::pair<std::size_t, std::size_t> CIPair;
    std::vector<CIPair> pairs;
    pairs.reserve(ls.size());

    // Create a sub dictionary containing all lines with location from ls.
    auto const &m = _dict.get<ByLocation>();
    for (auto l = ls.cbegin(); l != ls.cend(); ++l)
    {
        auto const p = m.equal_range(Line(*l));
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

    std::vector<std::size_t> global_indices;
    global_indices.reserve(pairs.size());
    for (auto p = pairs.cbegin(); p != pairs.cend(); ++p)
        global_indices.push_back(p->second);

    return global_indices;
}

}   // namespace AssemblerLib
