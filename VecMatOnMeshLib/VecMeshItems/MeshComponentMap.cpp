/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#include "MeshComponentMap.h"

#include <numeric>
#include <iostream>
#include <boost/range/algorithm/for_each.hpp>

#include "MeshLib/MeshSubsets.h"

namespace VecMatOnMeshLib
{

MeshComponentMap::MeshComponentMap(const std::vector<MeshLib::MeshSubsets*> &components, ComponentOrder order)
{
    // construct dict (and here we number global_index by component type)
    std::size_t global_index = 0;
    for (auto component = components.begin(); component != components.end(); ++component) {
        auto comp_id = std::distance(components.begin(), component);
        for (unsigned mesh_subset_index = 0; mesh_subset_index < (*component)->size(); mesh_subset_index++) {
			MeshLib::MeshSubset const& mesh_subset = (*component)->getMeshSubset(mesh_subset_index);
            std::size_t mesh_id = mesh_subset.getMeshID();
            // mesh items are ordered first by node, cell, ....
            for (std::size_t j=0; j<mesh_subset.getNNodes(); j++) {
                _dict.insert(MeshitemDataPosition(Location(mesh_id, MeshLib::MeshItemType::Node, j), comp_id, global_index++));
            }
            for (std::size_t j=0; j<mesh_subset.getNElements(); j++) {
                _dict.insert(MeshitemDataPosition(Location(mesh_id, MeshLib::MeshItemType::Cell, j), comp_id, global_index++));
            }
        }
    }

    if (order == ComponentOrder::BY_LOCATION)
        renumberByLocation();
}

void MeshComponentMap::renumberByLocation(std::size_t offset)
{
    std::size_t global_index = offset;

    auto &m = _dict.get<ByLocation>(); // view as sorted by mesh item
    for (auto itr_mesh_item=m.begin(); itr_mesh_item!=m.end(); ++itr_mesh_item) {
        MeshitemDataPosition pos = *itr_mesh_item;
        pos.global_index = global_index++;
        m.replace(itr_mesh_item, pos);
    }
}


std::size_t MeshComponentMap::getDataID(const Location &pos, unsigned compID) const
{
    auto &m = _dict.get<ByLocationAndComponent>();
    auto itr = m.find(MeshitemDataPosition(pos, compID, -1));
    return itr!=m.end() ? itr->global_index : -1;
}

std::vector<std::size_t> MeshComponentMap::getComponentIDs(const Location &pos) const
{
    auto &m = _dict.get<ByLocation>();
    auto p = m.equal_range(MeshitemDataPosition(pos, -1, -1));
    std::vector<std::size_t> vec_compID;
    for (auto itr=p.first; itr!=p.second; ++itr)
        vec_compID.push_back(itr->comp_id);
    return vec_compID;
}

std::vector<std::size_t> MeshComponentMap::getDataIDList(const Location &pos) const
{
    auto &m = _dict.get<ByLocation>();
    auto p = m.equal_range(MeshitemDataPosition(pos, -1, -1));
    std::vector<std::size_t> vec_dataID;
    for (auto itr=p.first; itr!=p.second; ++itr)
        vec_dataID.push_back(itr->global_index);
    return vec_dataID;
}


template <>
std::vector<std::size_t>
MeshComponentMap::getDataIDList<ComponentOrder::BY_LOCATION>(
    const std::vector<Location> &vec_pos) const
{
    // Create vector of global indices sorted by location containing all
    // locations given in vec_pos parameter.

    std::vector<std::size_t> vec_dataID;
    vec_dataID.reserve(vec_pos.size());

    auto &m = _dict.get<ByLocation>();
    for (auto &location : vec_pos) {
        auto p = m.equal_range(MeshitemDataPosition(location, -1, -1));
        for (auto itr=p.first; itr!=p.second; ++itr)
            vec_dataID.push_back(itr->global_index);
    }

    return vec_dataID;
}

template <>
std::vector<std::size_t>
MeshComponentMap::getDataIDList<ComponentOrder::BY_COMPONENT>(
    const std::vector<Location> &vec_pos) const
{
    // vector of (Component, global Index) pairs.
    typedef std::pair<std::size_t, std::size_t> CIPair;
    std::vector<CIPair> pairs;
    pairs.reserve(vec_pos.size());

    // Create a sub dictionary containing all lines with location from vec_pos.
    auto const &m = _dict.get<ByLocation>();
    for (auto const &location : vec_pos) {
        auto const p = m.equal_range(MeshitemDataPosition(location, -1, -1));
        for (auto itr=p.first; itr!=p.second; ++itr)
            pairs.emplace_back(itr->comp_id, itr->global_index);
    }

    auto CIPairLess = [](CIPair const& a, CIPair const& b)
        {
            return (a.first < b.first);
        };

    // Create vector of global indices from sub dictionary sorting by component.
    if (!std::is_sorted(pairs.begin(), pairs.end(), CIPairLess))
        std::stable_sort(pairs.begin(), pairs.end(), CIPairLess);

    std::vector<std::size_t> vec_dataID;
    vec_dataID.reserve(pairs.size());
    for (auto p = pairs.cbegin(); p != pairs.cend(); ++p)
        vec_dataID.push_back(p->second);

    return vec_dataID;
}

#ifndef NDEBUG
void disp(const MeshitemDataPosition &dat)
{
    dat.print();
}

void MeshComponentMap::print()
{

    boost::for_each(_dict, disp);
}
#endif	// NDEBUG

} // VecMatOnMeshLib

