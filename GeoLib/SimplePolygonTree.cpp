/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SimplePolygonTree.cpp
 *
 * Created on 2010-06-22 by Thomas Fischer
 */

#include "SimplePolygonTree.h"

namespace GeoLib {

SimplePolygonTree::SimplePolygonTree(const Polygon* polygon, SimplePolygonTree* parent) :
	_node (polygon), _parent (parent)
{}

SimplePolygonTree::~SimplePolygonTree()
{}

const Polygon* SimplePolygonTree::getPolygon () const
{
	return _node;
}

bool SimplePolygonTree::isPolygonInside (const SimplePolygonTree* polygon_hierarchy) const
{
	const Polygon* polygon (polygon_hierarchy->getPolygon());
	// check *all* points of polygon
	size_t n_pnts_polygon (polygon->getNumberOfPoints() - 1), cnt(0);
	for (size_t k(0); k<n_pnts_polygon && cnt == k; k++) {
		if (_node->isPntInPolygon (*(polygon->getPoint(k)))) {
			cnt++;
		}
	}
	// all points of the given polygon are contained in the
	if (cnt == n_pnts_polygon) return true;
	else {
		return false;
	}
}

void SimplePolygonTree::insertSimplePolygonTree (SimplePolygonTree* polygon_hierarchy)
{
	const Polygon* polygon (polygon_hierarchy->getPolygon());
	bool nfound (true);
	for (std::list<SimplePolygonTree*>::const_iterator it (_childs.begin());
		it != _childs.end() && nfound; it++)
	{
		// check all points of polygon
		size_t n_pnts_polygon (polygon->getNumberOfPoints()), cnt(0);
		for (size_t k(0); k<n_pnts_polygon && cnt == k; k++) {
			if (((*it)->getPolygon())->isPntInPolygon (*(polygon->getPoint(k))))
				cnt++;
		}
		// all points of the given polygon are contained in the current polygon
		if (cnt == n_pnts_polygon) {
			(*it)->insertSimplePolygonTree (polygon_hierarchy);
			nfound = false;
		}
	}
	if (nfound)
		_childs.push_back (polygon_hierarchy);
}

bool SimplePolygonTree::isGeoObjInside (const GeoObject* geo_obj) const
{
	if (dynamic_cast<const Point*>(geo_obj))
		return _node->isPntInPolygon (*(dynamic_cast<const Point*>(geo_obj)));

	if (dynamic_cast<const Polyline*>(geo_obj))
		return isPolylineInside (dynamic_cast<const Polyline*>(geo_obj));

	return false;
}

bool SimplePolygonTree::isPolylineInside (const Polyline* ply) const
{
	// check *all* points of polyline
	size_t n_pnts_polyline (ply->getNumberOfPoints() - 1), cnt(0);
	for (size_t k(0); k<n_pnts_polyline && cnt == k; k++) {
		if (_node->isPntInPolygon (*(ply->getPoint(k)))) {
			cnt++;
		}
	}
	// all points of the given polyline are contained in the polygon
	if (cnt == n_pnts_polyline) return true;

	return false;
}

void SimplePolygonTree::insertGeoObj (const GeoObject* geo_obj)
{
	// check if the geo object is contained in a child of this node
	bool nfound (true);
	for (std::list<SimplePolygonTree*>::const_iterator it (_childs.begin());
		it != _childs.end() && nfound; it++)
	{
		// check Point
		if (dynamic_cast<const Point*>(geo_obj)) {
			if (((*it)->getPolygon())->isPntInPolygon (*(dynamic_cast<const Point*>(geo_obj)))) {
				(*it)->insertGeoObj (geo_obj);
				nfound = false;
			}
		}
		// check Polyline
		if (nfound && dynamic_cast<const Polyline*>(geo_obj)) {
			const Polyline* ply (dynamic_cast<const Polyline*>(geo_obj));
			size_t n_pnts_polyline (ply->getNumberOfPoints()), cnt(0);
			// check all points of Polyline
			for (size_t k(0); k<n_pnts_polyline && cnt == k; k++) {
				if (((*it)->getPolygon())->isPntInPolygon (*(ply->getPoint(k))))
					cnt++;
			}
			// all points of the given polygon are contained in the current polygon
			if (cnt == n_pnts_polyline) {
				(*it)->insertGeoObj (geo_obj);
				nfound = false;
			}
		}
	}

	if (nfound) {
		_geo_objs.push_back (geo_obj);
	}

}

//void SimplePolygonTree::visitAndProcessNodes (FileIO::GMSHInterface& gmsh_io)
//{
//	if (getLevel() == 0) {
//		gmsh_io.writeGMSHPolygon (*_node);
//
//		std::list<SimplePolygonTree*>::iterator it (_childs.begin());
//		while (it != _childs.end()) {
//			(*it)->_visitAndProcessNodes (gmsh_io);
//			it++;
//		}
//		gmsh_io.writePlaneSurface ();
//	}
//}

//void SimplePolygonTree::_visitAndProcessNodes (FileIO::GMSHInterface& gmsh_io)
//{
//	gmsh_io.writeGMSHPolygon (*_node);
//
//	std::list<SimplePolygonTree*>::iterator it (_childs.begin());
//	while (it != _childs.end()) {
//		(*it)->_visitAndProcessNodes (gmsh_io);
//		it++;
//	}
//}

size_t SimplePolygonTree::getLevel () const
{
	if (_parent == NULL) return 0;
	else return 1+_parent->getLevel ();
}

void createPolygonTree (std::list<SimplePolygonTree*>& list_of_simple_polygon_hierarchies)
{
	std::list<SimplePolygonTree*>::iterator it0 (list_of_simple_polygon_hierarchies.begin()), it1;
	while (it0 != list_of_simple_polygon_hierarchies.end()) {
		it1 = it0;
		it1++;
		while (it1 != list_of_simple_polygon_hierarchies.end()) {
			if ((*it0)->isPolygonInside (*it1)) {
				(*it0)->insertSimplePolygonTree (*it1);
				it1 = list_of_simple_polygon_hierarchies.erase (it1);
			} else {
				if ((*it1)->isPolygonInside (*it0)) {
					(*it1)->insertSimplePolygonTree (*it0);
					(*it1)->insertSimplePolygonTree (*it0);
					it0 = list_of_simple_polygon_hierarchies.erase (it0);
				}

				it1++;
			}
		}
		it0++;
	}
}


} // end namespace GeoLib
