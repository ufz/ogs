/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-22
 * \brief  Implementation of the SimplePolygonTree class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SimplePolygonTree.h"

namespace GeoLib
{
SimplePolygonTree::SimplePolygonTree(Polygon* polygon,
                                     SimplePolygonTree* parent)
    : _node_polygon(polygon), _parent(parent)
{
}

SimplePolygonTree::~SimplePolygonTree()
{
    for (auto* child : _children)
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
    Polygon const& polygon = polygon_hierarchy->polygon();
    bool nfound(true);
    for (std::list<SimplePolygonTree*>::const_iterator it(_children.begin());
         it != _children.end() && nfound;
         ++it)
    {
        if ((*it)->polygon().isPolylineInPolygon(polygon))
        {
            (*it)->insertSimplePolygonTree(polygon_hierarchy);
            nfound = false;
        }
    }
    if (nfound)
    {
        _children.push_back(polygon_hierarchy);
        polygon_hierarchy->setParent(this);
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
