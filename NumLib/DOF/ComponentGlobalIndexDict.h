/**
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <limits>

#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index_container.hpp>
#include <utility>

#include "MeshLib/Location.h"
#include "NumLib/NumericsConfig.h"

namespace NumLib
{

/// \internal
namespace detail
{

struct Line
{
    MeshLib::Location location;

    // Physical component
    int comp_id;

    // Position in global matrix or vector
    GlobalIndexType global_index;

    Line(MeshLib::Location l, int c, GlobalIndexType i)
        : location(std::move(l)), comp_id(c), global_index(i)
    {}

    Line(MeshLib::Location l, int c)
        : location(std::move(l)),
          comp_id(c),
          global_index(std::numeric_limits<GlobalIndexType>::max())
    {}

    explicit Line(MeshLib::Location l)
        : location(std::move(l)),
          comp_id(std::numeric_limits<int>::max()),
          global_index(std::numeric_limits<GlobalIndexType>::max())
    {}

    friend std::ostream& operator<<(std::ostream& os, Line const& l)
    {
        return os << l.location << ", " << l.comp_id << ", " << l.global_index;
    }
};

struct LineByLocationComparator
{
    bool operator()(Line const& a, Line const& b) const
    {
        return a.location < b.location;
    }
};

struct LineByLocationAndComponentComparator
{
    bool operator()(Line const& a, Line const& b) const
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

using ComponentGlobalIndexDict = boost::multi_index::multi_index_container<
    Line, boost::multi_index::indexed_by<
              boost::multi_index::ordered_unique<
                  boost::multi_index::tag<ByLocationAndComponent>,
                  boost::multi_index::identity<Line>,
                  LineByLocationAndComponentComparator>,
              boost::multi_index::ordered_non_unique<
                  boost::multi_index::tag<ByLocation>,
                  boost::multi_index::identity<Line>, LineByLocationComparator>,
              boost::multi_index::ordered_non_unique<
                  boost::multi_index::tag<ByComponent>,
                  boost::multi_index::member<Line, int, &Line::comp_id>>,
              boost::multi_index::ordered_non_unique<
                  boost::multi_index::tag<ByGlobalIndex>,
                  boost::multi_index::member<Line, GlobalIndexType,
                                             &Line::global_index>>>>;

}    // namespace detail
}    // namespace NumLib
