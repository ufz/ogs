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


#ifndef MESHITEMDATAPOSITIONDICTIONARY_H_
#define MESHITEMDATAPOSITIONDICTIONARY_H_

#include <iostream>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/composite_key.hpp>

namespace VecMatOnMeshLib
{

struct MeshitemDataPosition
{
    Location location;

    // Physical component
    std::size_t comp_id;

    // Position in global matrix or vector
    std::size_t global_index;

    MeshitemDataPosition(Location const& location,
                         std::size_t comp_id,
                         std::size_t global_index)
    : location(location),
      comp_id(comp_id),
      global_index(global_index)
    {}

    void print() const
    {
        std::cout << location << ", " << comp_id << ", " << global_index << "\n";
    }
};

struct MeshitemDataPositionByLocationComparator
{
    bool operator()(MeshitemDataPosition const& a, MeshitemDataPosition const& b) const
    {
        return a.location < b.location;
    }
};

struct MeshitemDataPositionByLocationAndComponentComparator
{
    bool operator()(MeshitemDataPosition const& a, MeshitemDataPosition const& b) const
    {
        if (a.location < b.location)
            return true;
        if (b.location < a.location)
            return false;

        // a.loc == b.loc
        return a.comp_id < b.comp_id;
    }
};

struct ByLocation {};
struct ByLocationAndComponent {};
struct ByComponent {};
struct ByGlobalIndex {};

typedef boost::multi_index::multi_index_container<
        MeshitemDataPosition,
        boost::multi_index::indexed_by
        <
            boost::multi_index::ordered_unique
            <
                boost::multi_index::tag<ByLocationAndComponent>,
                boost::multi_index::identity<MeshitemDataPosition>,
                MeshitemDataPositionByLocationAndComponentComparator
            >,
            boost::multi_index::ordered_non_unique
            <
                boost::multi_index::tag<ByLocation>,
                boost::multi_index::identity<MeshitemDataPosition>,
                MeshitemDataPositionByLocationComparator
            >,
            boost::multi_index::ordered_non_unique
            <
                boost::multi_index::tag<ByComponent>,
                boost::multi_index::member<MeshitemDataPosition, std::size_t, &MeshitemDataPosition::comp_id>
            >,
            boost::multi_index::ordered_non_unique
            <
                boost::multi_index::tag<ByGlobalIndex>,
                boost::multi_index::member<MeshitemDataPosition, std::size_t, &MeshitemDataPosition::global_index>
            >
        >
    > MeshitemDataPositionDictionary;

} // VecMatOnMeshLib

#endif /* MESHITEMDATAPOSITIONDICTIONARY_H_ */
