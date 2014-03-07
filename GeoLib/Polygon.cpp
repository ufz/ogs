/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-21
 * \brief  Implementation of the Polygon class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cmath>
#include <cstdlib> // for exit

// ThirdParty/logog
#include "logog/include/logog.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "Polygon.h"

// MathLib
#include "AnalyticalGeometry.h"
#include "MathTools.h"
#include "Vector3.h"

// BaseLib
#include "quicksort.h"

namespace GeoLib
{
Polygon::Polygon(const Polyline &ply, bool init) :
	Polyline(ply), _aabb(ply.getPointsVec(), ply._ply_pnt_ids)
{
	if (init)
		initialise ();
}

Polygon::~Polygon()
{
	// remove polygons from list
	for (std::list<Polygon*>::iterator it (_simple_polygon_list.begin());
	     it != _simple_polygon_list.end(); it++)
		// the first entry of the list can be a pointer the object itself
		if (*it != this)
			delete *it;
}

bool Polygon::initialise ()
{
	if (this->isClosed()) {
		ensureCWOrientation();
		return true;
	} else {
		WARN("Polygon::initialise(): base polyline is not closed.");
		return false;
	}
}

bool Polygon::isPntInPolygon (GeoLib::Point const & pnt) const
{
	GeoLib::Point const& min_aabb_pnt (_aabb.getMinPoint());
	GeoLib::Point const& max_aabb_pnt (_aabb.getMaxPoint());

	if (pnt[0] < min_aabb_pnt[0] || max_aabb_pnt[0] < pnt[0] || pnt[1] < min_aabb_pnt[1] ||
	    max_aabb_pnt[1] < pnt[1])
		return false;

	std::size_t n_intersections (0);
	GeoLib::Point s;

	if (_simple_polygon_list.empty ()) {
		const std::size_t n_nodes (getNumberOfPoints() - 1);
		for (std::size_t k(0); k < n_nodes; k++) {
			if (((*(getPoint(k)))[1] <= pnt[1] && pnt[1] <= (*(getPoint(k + 1)))[1]) ||
			    ((*(getPoint(k + 1)))[1] <= pnt[1] && pnt[1] <= (*(getPoint(k)))[1])) {
				switch (getEdgeType(k, pnt))
				{
				case EdgeType::TOUCHING:
					return true;
				case EdgeType::CROSSING:
					n_intersections++;
					break;
				case EdgeType::INESSENTIAL:
					break;
				default:
					// do nothing
					;
				}
			}
		}
		if (n_intersections % 2 == 1)
			return true;
	} else {
		for (std::list<Polygon*>::const_iterator it (_simple_polygon_list.begin());
		     it != _simple_polygon_list.end(); ++it) {
			if ((*it)->isPntInPolygon (pnt))
				return true;
		}
		return false;
	}
	return false;
}

bool Polygon::isPntInPolygon(double x, double y, double z) const
{
	const GeoLib::Point pnt(x,y,z);
	return isPntInPolygon (pnt);
}

bool Polygon::isPolylineInPolygon(const Polyline& ply) const
{
	std::size_t ply_size (ply.getNumberOfPoints()), cnt (0);
	for (std::size_t k(0); k < ply_size; k++) {
		if (isPntInPolygon (*(ply.getPoint(k)))) {
			cnt++;
		}
	}

	if (cnt == ply_size)
		return true;
	return false;
}

bool Polygon::isPartOfPolylineInPolygon(const Polyline& ply) const
{
	const std::size_t ply_size (ply.getNumberOfPoints());
	// check points
	for (std::size_t k(0); k < ply_size; k++) {
		if (isPntInPolygon (*(ply.getPoint(k)))) {
			return true;
		}
	}
	// check segment intersections
	GeoLib::Point* s (new GeoLib::Point (0,0,0));
	const std::size_t n_nodes(getNumberOfPoints() - 1);
	for (std::size_t k(0); k < ply_size - 1; k++) {
		for (std::size_t j(0); j < n_nodes; j++) {
			if (GeoLib::lineSegmentIntersect(*(getPoint(j)), *(getPoint(j + 1)),
							*(ply.getPoint(k)), *(ply.getPoint(k + 1)), *s)) {
				delete s;
				return true;
			}
		}
	}

	delete s;
	return false;
}

bool Polygon::getNextIntersectionPointPolygonLine (GeoLib::Point const & a,
                GeoLib::Point const & b, GeoLib::Point* intersection_pnt,
                std::size_t& seg_num) const
{
	if (_simple_polygon_list.empty()) {
		const std::size_t n_segments(getNumberOfPoints() - 1);
		for (std::size_t k(seg_num); k < n_segments; k++) {
			if (GeoLib::lineSegmentIntersect(*(getPoint(k)), *(getPoint(k + 1)), a, b, *intersection_pnt)) {
				seg_num = k;
				return true;
			}
		}
	} else {
		for (std::list<Polygon*>::const_iterator it(_simple_polygon_list.begin()); it
					!= _simple_polygon_list.end(); ++it) {
			const Polygon* polygon(*it);
			const std::size_t n_nodes_simple_polygon(polygon->getNumberOfPoints() - 1);
			for (std::size_t k(0); k < n_nodes_simple_polygon; k++) {
				if (GeoLib::lineSegmentIntersect(*(polygon->getPoint(k)), *(polygon->getPoint(k + 1)),
								a, b, *intersection_pnt)) {
					seg_num = k;
					return true;
				}
			}
		}
	}
	return false;
}

const std::list<Polygon*>& Polygon::getListOfSimplePolygons()
{
	if (_simple_polygon_list.empty())
		_simple_polygon_list.push_back (this);
	return _simple_polygon_list;
}

void Polygon::computeListOfSimplePolygons ()
{
	if (!_simple_polygon_list.empty())
		return;

	_simple_polygon_list.push_back (this);
	splitPolygonAtPoint (_simple_polygon_list.begin());
	splitPolygonAtIntersection (_simple_polygon_list.begin());

	for (std::list<Polygon*>::iterator it (_simple_polygon_list.begin());
	     it != _simple_polygon_list.end(); it++)
		(*it)->initialise ();
}

EdgeType Polygon::getEdgeType (std::size_t k, GeoLib::Point const & pnt) const
{
	switch (getLocationOfPoint(k, pnt))
	{
	case Location::LEFT: {
		const GeoLib::Point & v (*(getPoint(k)));
		const GeoLib::Point & w (*(getPoint(k + 1)));
		if (v[1] < pnt[1] && pnt[1] <= w[1])
			return EdgeType::CROSSING;
		else
			return EdgeType::INESSENTIAL;
	}
	case Location::RIGHT: {
		const GeoLib::Point & v (*(getPoint(k)));
		const GeoLib::Point & w (*(getPoint(k + 1)));
		if (w[1] < pnt[1] && pnt[1] <= v[1])
			return EdgeType::CROSSING;
		else
			return EdgeType::INESSENTIAL;
	}
	case Location::BETWEEN:
	case Location::SOURCE:
	case Location::DESTINATION:
		return EdgeType::TOUCHING;
	default:
		return EdgeType::INESSENTIAL;
	}
}

void Polygon::ensureCWOrientation ()
{
	// *** pre processing: rotate points to xy-plan
	// *** copy points to vector - last point is identical to the first
	std::size_t n_pnts (this->getNumberOfPoints() - 1);
	std::vector<GeoLib::Point*> tmp_polygon_pnts;
	for (std::size_t k(0); k < n_pnts; k++)
		tmp_polygon_pnts.push_back (new GeoLib::Point (*(this->getPoint(k))));

	// *** calculate supporting plane (plane normal and
	MathLib::Vector3 plane_normal;
	double d;
	GeoLib::getNewellPlane(tmp_polygon_pnts, plane_normal, d);

	// *** rotate if necessary
	double tol (sqrt(std::numeric_limits<double>::min()));
	if (fabs(plane_normal[0]) > tol || fabs(plane_normal[1]) > tol)
		// rotate copied points into x-y-plane
		GeoLib::rotatePointsToXY(plane_normal, tmp_polygon_pnts);

	for (std::size_t k(0); k < tmp_polygon_pnts.size(); k++)
		(*(tmp_polygon_pnts[k]))[2] = 0.0; // should be -= d but there are numerical errors

	// *** get the left most upper point
	std::size_t min_x_max_y_idx (0); // for orientation check
	for (std::size_t k(0); k < n_pnts; k++)
		if ((*(tmp_polygon_pnts[k]))[0] <= (*(tmp_polygon_pnts[min_x_max_y_idx]))[0])
		{
			if ((*(tmp_polygon_pnts[k]))[0] < (*(tmp_polygon_pnts[min_x_max_y_idx]))[0])
				min_x_max_y_idx = k;
			else if ((*(tmp_polygon_pnts[k]))[1] >
			         (*(tmp_polygon_pnts[min_x_max_y_idx]))[1])
				min_x_max_y_idx = k;

		}
	// *** determine orientation
	GeoLib::Orientation orient;
	if (0 < min_x_max_y_idx && min_x_max_y_idx < n_pnts - 2)
		orient = GeoLib::getOrientation (
		        tmp_polygon_pnts[min_x_max_y_idx - 1],
		        tmp_polygon_pnts[min_x_max_y_idx],
		        tmp_polygon_pnts[min_x_max_y_idx + 1]);
	else
	{
		if (0 == min_x_max_y_idx)
			orient = GeoLib::getOrientation (
			        tmp_polygon_pnts[n_pnts - 1],
			        tmp_polygon_pnts[0],
			        tmp_polygon_pnts[1]);
		else
			orient = GeoLib::getOrientation (
			        tmp_polygon_pnts[n_pnts - 2],
			        tmp_polygon_pnts[n_pnts - 1],
			        tmp_polygon_pnts[0]);
	}

	if (orient == GeoLib::CCW)
	{
		// switch orientation
		std::size_t tmp_n_pnts (n_pnts);
		tmp_n_pnts++; // include last point of polygon (which is identical to the first)
		for (std::size_t k(0); k < tmp_n_pnts / 2; k++)
			std::swap (_ply_pnt_ids[k], _ply_pnt_ids[tmp_n_pnts - 1 - k]);
	}

	for (std::size_t k(0); k < n_pnts; k++)
		delete tmp_polygon_pnts[k];
}

void Polygon::splitPolygonAtIntersection (std::list<Polygon*>::iterator polygon_it)
{
	std::size_t idx0 (0), idx1 (0);
	while (polygon_it != _simple_polygon_list.end())
	{
		GeoLib::Point* intersection_pnt (new GeoLib::Point);
		bool is_simple (!GeoLib::lineSegmentsIntersect (*polygon_it,
		                                                 idx0,
		                                                 idx1,
		                                                 *intersection_pnt));
		if (!is_simple)
		{
			// adding intersection point to pnt_vec
			std::size_t intersection_pnt_id (_ply_pnts.size());
			const_cast<std::vector<Point*>& >(_ply_pnts).push_back (intersection_pnt);

			// split Polygon
			if (idx0 > idx1)
				std::swap (idx0, idx1);

			GeoLib::Polygon* polygon0 (new GeoLib::Polygon(
			                                   (*polygon_it)->getPointsVec(), false));
			for (std::size_t k(0); k <= idx0; k++)
				polygon0->addPoint ((*polygon_it)->getPointID (k));
			polygon0->addPoint (intersection_pnt_id);
			for (std::size_t k(idx1 + 1); k < (*polygon_it)->getNumberOfPoints(); k++)
				polygon0->addPoint ((*polygon_it)->getPointID (k));
			if (!polygon0->initialise())
			{
				ERR("Polygon::splitPolygonAtIntersection(): Initialization of polygon0 failed.");
				exit (1);
			}

			GeoLib::Polygon* polygon1 (new GeoLib::Polygon(
			                                   (*polygon_it)->getPointsVec(), false));
			polygon1->addPoint (intersection_pnt_id);
			for (std::size_t k(idx0 + 1); k <= idx1; k++)
				polygon1->addPoint ((*polygon_it)->getPointID (k));
			polygon1->addPoint (intersection_pnt_id);
			if (!polygon1->initialise())
			{
				ERR("Polygon::splitPolygonAtIntersection(): Initialization of polygon1 failed.");
				exit (1);
			}

			// remove original polyline and add two new polylines
			std::list<GeoLib::Polygon*>::iterator polygon0_it, polygon1_it;
			polygon_it = _simple_polygon_list.erase (polygon_it);
			polygon1_it = _simple_polygon_list.insert (polygon_it, polygon1);
			polygon0_it = _simple_polygon_list.insert (polygon1_it, polygon0);

			splitPolygonAtIntersection (polygon0_it);
			splitPolygonAtIntersection (polygon1_it);
		}
		else
			delete intersection_pnt;
		++polygon_it;
	}
}

void Polygon::splitPolygonAtPoint (std::list<GeoLib::Polygon*>::iterator polygon_it)
{
	std::size_t n ((*polygon_it)->getNumberOfPoints() - 1), idx0 (0), idx1(0);
	std::size_t* id_vec (new std::size_t[n]), *perm (new std::size_t[n]);
	for (std::size_t k(0); k < n; k++)
	{
		id_vec[k] = (*polygon_it)->getPointID (k);
		perm[k] = k;
	}

	BaseLib::quicksort (id_vec, 0, n, perm);

	for (std::size_t k(0); k < n - 1; k++)
		if (id_vec[k] == id_vec[k + 1])
		{
			idx0 = perm[k];
			idx1 = perm[k + 1];
			delete [] perm;
			delete [] id_vec;

			if (idx0 > idx1)
				std::swap (idx0, idx1);

			// create two closed polylines
			GeoLib::Polygon* polygon0 (new GeoLib::Polygon(*(*polygon_it)));
			for (std::size_t k(0); k <= idx0; k++)
				polygon0->addPoint ((*polygon_it)->getPointID (k));
			for (std::size_t k(idx1 + 1); k < (*polygon_it)->getNumberOfPoints(); k++)
				polygon0->addPoint ((*polygon_it)->getPointID (k));
			polygon0->initialise();

			GeoLib::Polygon* polygon1 (new GeoLib::Polygon(*(*polygon_it)));
			for (std::size_t k(idx0); k <= idx1; k++)
				polygon1->addPoint ((*polygon_it)->getPointID (k));
			polygon1->initialise();

			// remove original polygon and add two new polygons
			std::list<GeoLib::Polygon*>::iterator polygon0_it, polygon1_it;
			polygon1_it =
			        _simple_polygon_list.insert (_simple_polygon_list.erase (
			                                             polygon_it), polygon1);
			polygon0_it = _simple_polygon_list.insert (polygon1_it, polygon0);

			splitPolygonAtPoint (polygon0_it);
			splitPolygonAtPoint (polygon1_it);

			return;
		}
	delete [] perm;
	delete [] id_vec;
}

GeoLib::Polygon* createPolygonFromCircle (GeoLib::Point const& middle_pnt, double radius,
                                          std::vector<GeoLib::Point*> & pnts, std::size_t resolution)
{
	const std::size_t off_set (pnts.size());
	// create points
	double angle (2.0 * M_PI / resolution);
	for (std::size_t k(0); k < resolution; k++)
	{
		GeoLib::Point* pnt (new GeoLib::Point(middle_pnt.getCoords()));
		(*pnt)[0] += radius * cos (k * angle);
		(*pnt)[1] += radius * sin (k * angle);
		pnts.push_back (pnt);
	}

	// create polygon
	GeoLib::Polygon* polygon (new GeoLib::Polygon (pnts, false));
	for (std::size_t k(0); k < resolution; k++)
		polygon->addPoint (k + off_set);
	polygon->addPoint (off_set);

	return polygon;
}

bool operator==(Polygon const& lhs, Polygon const& rhs)
{
	if (lhs.getNumberOfPoints() != rhs.getNumberOfPoints())
		return false;

	const std::size_t n(lhs.getNumberOfPoints());
	const std::size_t start_pnt(lhs.getPointID(0));

	// search start point of first polygon in second polygon
	bool nfound(true);
	std::size_t k(0);
	for (; k < n-1 && nfound; k++) {
		if (start_pnt == rhs.getPointID(k)) {
			nfound = false;
			break;
		}
	}

	// case: start point not found in second polygon
	if (nfound) return false;

	// *** determine direction
	// opposite direction
	if (k == n-2) {
		for (k=1; k<n-1; k++) {
			if (lhs.getPointID(k) != rhs.getPointID(n-1-k)) {
				return false;
			}
		}
		return true;
	}

	// same direction - start point of first polygon at arbitrary position in second polygon
	if (lhs.getPointID(1) == rhs.getPointID(k+1)) {
		std::size_t j(k+2);
		for (; j<n-1; j++) {
			if (lhs.getPointID(j-k) != rhs.getPointID(j)) {
				return false;
			}
		}
		j=0; // new start point at second polygon
		for (; j<k+1; j++) {
			if (lhs.getPointID(n-(k+2)+j+1) != rhs.getPointID(j)) {
				return false;
			}
		}
		return true;
	} else {
		// opposite direction with start point of first polygon at arbitrary position
		// *** ATTENTION
		WARN("operator==(Polygon const& lhs, Polygon const& rhs) - not tested case (implementation is probably buggy) - please contact thomas.fischer@ufz.de mentioning the problem.");
		// in second polygon
		if (lhs.getPointID(1) == rhs.getPointID(k-1)) {
			std::size_t j(k-2);
			for (; j>0; j--) {
				if (lhs.getPointID(k-2-j) != rhs.getPointID(j)) {
					return false;
				}
			}
			// new start point at second polygon - the point n-1 of a polygon is equal to the
			// first point of the polygon (for this reason: n-2)
			j=n-2;
			for (; j>k-1; j--) {
				if (lhs.getPointID(n-2+j+k-2) != rhs.getPointID(j)) {
					return false;
				}
			}
			return true;
		} else {
			return false;
		}
	}
}

} // end namespace GeoLib
