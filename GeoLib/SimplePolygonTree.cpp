/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-22
 * \brief  Implementation of the SimplePolygonTree class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SimplePolygonTree.h"

namespace GeoLib
{
SimplePolygonTree::SimplePolygonTree(Polygon * polygon, SimplePolygonTree * parent) :
    node_polygon_ (polygon), parent_ (parent)
{}

SimplePolygonTree::~SimplePolygonTree()
{
    for (auto * child : children_) {
        delete child;
    }
}

bool SimplePolygonTree::isPolygonInside (const SimplePolygonTree* polygon_hierarchy) const
{
    return node_polygon_->isPolylineInPolygon(*(polygon_hierarchy->getPolygon()));
}

void SimplePolygonTree::insertSimplePolygonTree (SimplePolygonTree* polygon_hierarchy)
{
    const Polygon* polygon (polygon_hierarchy->getPolygon());
    bool nfound (true);
    for (std::list<SimplePolygonTree*>::const_iterator it (children_.begin());
         it != children_.end() && nfound; ++it) {
        if (((*it)->getPolygon())->isPolylineInPolygon (*(polygon))) {
            (*it)->insertSimplePolygonTree (polygon_hierarchy);
            nfound = false;
        }
    }
    if (nfound) {
        children_.push_back (polygon_hierarchy);
        polygon_hierarchy->setParent(this);
    }
}

const Polygon* SimplePolygonTree::getPolygon () const
{
    return node_polygon_;
}

} // end namespace GeoLib
