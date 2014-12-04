/**
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <iostream>

#include "MeshLib/MeshSubsets.h"

#ifdef USE_PETSC
#include "MeshLib/NodePartitionedMesh.h"
#endif

#include "MeshComponentMap.h"

namespace AssemblerLib
{

using namespace detail;

template <typename T_INT_TYPE>
std::vector<T_INT_TYPE>
MeshComponentMap::getGlobalIndicesByLocation(
    std::vector<Location> const &ls) const
{
    // Create vector of global indices sorted by location containing all
    // locations given in ls parameter.

    std::vector<T_INT_TYPE> global_indices;
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

template <typename T_INT_TYPE>
std::vector<T_INT_TYPE>
MeshComponentMap::getGlobalIndicesByComponent(
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

    std::vector<T_INT_TYPE> global_indices;
    global_indices.reserve(pairs.size());
    for (auto p = pairs.cbegin(); p != pairs.cend(); ++p)
        global_indices.push_back(p->second);

    return global_indices;
}

}   // namespace AssemblerLib
