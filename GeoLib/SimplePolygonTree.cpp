/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-22
 * \brief  Implementation of the SimplePolygonTree class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SimplePolygonTree.h"

namespace GeoLib
{
SimplePolygonTree::SimplePolygonTree(Polygon * polygon, SimplePolygonTree * parent) :
    _node_polygon (polygon), _parent (parent)
{}

SimplePolygonTree::~SimplePolygonTree()
{
    for (auto * child : _children) {
        delete child;
    }
}

bool SimplePolygonTree::isPolygonInside (const SimplePolygonTree* polygon_hierarchy) const
{
    return _node_polygon->isPolylineInPolygon(*(polygon_hierarchy->getPolygon()));
}

void SimplePolygonTree::insertSimplePolygonTree (SimplePolygonTree* polygon_hierarchy)
{
    const Polygon* polygon (polygon_hierarchy->getPolygon());
    bool nfound (true);
    for (std::list<SimplePolygonTree*>::const_iterator it (_children.begin());
         it != _children.end() && nfound; ++it) {
        if (((*it)->getPolygon())->isPolylineInPolygon (*(polygon))) {
            (*it)->insertSimplePolygonTree (polygon_hierarchy);
            nfound = false;
        }
    }
    if (nfound) {
        _children.push_back (polygon_hierarchy);
        polygon_hierarchy->setParent(this);
    }
}

const Polygon* SimplePolygonTree::getPolygon () const
{
    return _node_polygon;
}

} // end namespace GeoLib
