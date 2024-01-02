/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-22
 * \brief  Implementation of the SimplePolygonTree class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SimplePolygonTree.h"

#include <range/v3/algorithm/find_if.hpp>

namespace GeoLib
{
SimplePolygonTree::SimplePolygonTree(Polygon* polygon,
                                     SimplePolygonTree* parent)
    : _node_polygon(polygon), _parent(parent)
{
}

SimplePolygonTree::~SimplePolygonTree()
{
    for (auto const* child : _children)
    {
        delete child;
    }
}

bool SimplePolygonTree::isRoot() const
{
    return _parent == nullptr;
}

bool SimplePolygonTree::isPolygonInside(
    const SimplePolygonTree* polygon_hierarchy) const
{
    return _node_polygon->isPolylineInPolygon(polygon_hierarchy->polygon());
}

const SimplePolygonTree* SimplePolygonTree::parent() const
{
    return _parent;
}

void SimplePolygonTree::insertSimplePolygonTree(
    SimplePolygonTree* polygon_hierarchy)
{
    auto const child = ranges::find_if(
        _children, [&p = polygon_hierarchy->polygon()](auto const* c)
        { return c->polygon().isPolylineInPolygon(p); });

    if (child != std::end(_children))
    {
        (*child)->insertSimplePolygonTree(polygon_hierarchy);
    }
    else
    {
        _children.push_back(polygon_hierarchy);
        polygon_hierarchy->_parent = this;
    }
}

Polygon& SimplePolygonTree::polygon()
{
    return *_node_polygon;
}
Polygon const& SimplePolygonTree::polygon() const
{
    return *_node_polygon;
}

}  // end namespace GeoLib
