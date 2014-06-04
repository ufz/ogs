/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-22
 * \brief  Implementation of the SimplePolygonTree class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
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
	for (std::list<SimplePolygonTree*>::const_iterator it (_childs.begin());
		     it != _childs.end(); ++it) {
		delete *it;
	}
}

bool SimplePolygonTree::isPolygonInside (const SimplePolygonTree* polygon_hierarchy) const
{
	const Polygon* polygon (polygon_hierarchy->getPolygon());
	// check *all* points of polygon
	size_t n_pnts_polygon (polygon->getNumberOfPoints() - 1), cnt(0);
	for (size_t k(0); k < n_pnts_polygon && cnt == k; k++) {
		if (_node_polygon->isPntInPolygon (*(polygon->getPoint(k)))) {
			cnt++;
		}
	}

	// all points of the given polygon are contained in the
	if (cnt == n_pnts_polygon)
		return true;
	else
		return false;
}

void SimplePolygonTree::insertSimplePolygonTree (SimplePolygonTree* polygon_hierarchy)
{
	const Polygon* polygon (polygon_hierarchy->getPolygon());
	bool nfound (true);
	for (std::list<SimplePolygonTree*>::const_iterator it (_childs.begin());
	     it != _childs.end() && nfound; ++it) {
		// check all points of polygon
		size_t n_pnts_polygon (polygon->getNumberOfPoints()), cnt(0);
		for (size_t k(0); k < n_pnts_polygon && cnt == k; k++) {
			if (((*it)->getPolygon())->isPntInPolygon (*(polygon->getPoint(k)))) {
				cnt++;
			}
		}
		// all points of the given polygon are contained in the current polygon
		if (cnt == n_pnts_polygon) {
			(*it)->insertSimplePolygonTree (polygon_hierarchy);
			nfound = false;
		}
	}
	if (nfound) {
		_childs.push_back (polygon_hierarchy);
		polygon_hierarchy->setParent(this);
	}
}

const Polygon* SimplePolygonTree::getPolygon () const
{
	return _node_polygon;
}

} // end namespace GeoLib
