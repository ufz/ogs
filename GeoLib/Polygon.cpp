/*
 * Polygon.cpp
 *
 *  Created on: Jun 21, 2010
 *      Author: TF
 */

#include <cstdlib> // for exit
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "Polygon.h"

// MathLib
#include "AnalyticalGeometry.h"
#include "MathTools.h"
#include "Vector3.h"

// Base
#include "quicksort.h"
#include "swap.h"

namespace GEOLIB {

Polygon::Polygon(const Polyline &ply, bool init) :
	Polyline (ply)
{
	if (init) {
		initialise ();
	}
}

Polygon::Polygon (const std::vector<Point*>& pnt_vec) :
	Polyline (pnt_vec)
{}

Polygon::~Polygon()
{
	// remove polygons from list
	for (std::list<Polygon*>::iterator it (_simple_polygon_list.begin()); it != _simple_polygon_list.end(); it++) {
		// the first entry of the list can be a pointer the object itself
		if (*it != this)
			delete *it;
	}
}

bool Polygon::initialise ()
{
	if (this->isClosed()) {
		calculateAxisAlignedBoundingBox();
		ensureCWOrientation();
		return true;
	} else {
		std::cerr << "ERROR in Polygon::initialise() - base polyline is not closed" << std::endl;
		return false;
	}
}

bool Polygon::isPntInPolygon (GEOLIB::Point const & pnt) const
{
	GEOLIB::Point min_aabb_pnt (_aabb.getMinPoint());
	GEOLIB::Point max_aabb_pnt (_aabb.getMaxPoint());

	if (pnt[0] < min_aabb_pnt[0] || max_aabb_pnt[0] < pnt[0] || pnt[1] < min_aabb_pnt[1] || max_aabb_pnt[1] < pnt[1])
		return false;

	size_t n_intersections (0);
	GEOLIB::Point s;

	if (_simple_polygon_list.empty ()) {
		const size_t n_nodes (getNumberOfPoints()-1);
		for (size_t k(0); k<n_nodes; k++) {
			if (((*(getPoint(k)))[1] <= pnt[1] && pnt[1] <= (*(getPoint(k+1)))[1]) ||
					((*(getPoint(k+1)))[1] <= pnt[1] && pnt[1] <= (*(getPoint(k)))[1])) {
				switch (getEdgeType(k, pnt)) {
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
		if (n_intersections%2 == 1) return true;
	} else {
		for (std::list<Polygon*>::const_iterator it (_simple_polygon_list.begin());
			it != _simple_polygon_list.end(); ++it) {
			if ((*it)->isPntInPolygon (pnt)) return true;
		}
		return false;
	}
	return false;
}

bool Polygon::isPntInPolygon (double x, double y, double z) const
{
	const GEOLIB::Point pnt(x,y,z);
	return isPntInPolygon (pnt);
}

bool Polygon::isPolylineInPolygon (const Polyline& ply) const
{
	size_t ply_size (ply.getNumberOfPoints()), cnt (0);
	for (size_t k(0); k<ply_size; k++) {
		if (isPntInPolygon (*(ply[k]))) {
			cnt++;
		}
	}
	if (cnt == ply_size)
		return true;
	return false;
}

GEOLIB::Point* Polygon::getIntersectionPointPolygonLine (GEOLIB::Point const & a, GEOLIB::Point const & b) const
{
	GEOLIB::Point* s (new GEOLIB::Point (0,0,0));

	if (_simple_polygon_list.empty ()) {
		const size_t n_nodes (getNumberOfPoints()-1);
		for (size_t k(0); k<n_nodes; k++) {
			if (MathLib::lineSegmentIntersect (*(getPoint(k)), *(getPoint(k+1)), a, b, *s)) {
				return s;
			}
		}
	} else {
		for (std::list<Polygon*>::const_iterator it (_simple_polygon_list.begin());
			it != _simple_polygon_list.end(); ++it) {
			const Polygon* polygon (*it);
			const size_t n_nodes_simple_polygon (polygon->getNumberOfPoints()-1);
			for (size_t k(0); k<n_nodes_simple_polygon; k++) {
				if (MathLib::lineSegmentIntersect (*(polygon->getPoint(k)), *(polygon->getPoint(k+1)), a, b, *s)) {
					return s;
				}
			}
		}
	}
	delete s;
	return NULL;
}

const std::list<Polygon*>& Polygon::getListOfSimplePolygons()
{
	if (_simple_polygon_list.empty())
		_simple_polygon_list.push_back (this);
	return _simple_polygon_list;
}

void Polygon::computeListOfSimplePolygons ()
{
	if (! _simple_polygon_list.empty())
		return;

	_simple_polygon_list.push_back (this);
	splitPolygonAtPoint (_simple_polygon_list.begin());
	splitPolygonAtIntersection (_simple_polygon_list.begin());

	for (std::list<Polygon*>::iterator it (_simple_polygon_list.begin());
		it != _simple_polygon_list.end(); it++) {
		(*it)->initialise ();
	}
}

EdgeType::value Polygon::getEdgeType (size_t k, GEOLIB::Point const & pnt) const
{
	switch (getLocationOfPoint(k, pnt)) {
	case Location::LEFT: {
		const GEOLIB::Point & v (*(getPoint(k)));
		const GEOLIB::Point & w (*(getPoint(k+1)));
		if (v[1] < pnt[1] && pnt[1] <= w[1]) return EdgeType::CROSSING;
		else return EdgeType::INESSENTIAL;
		break;
	}
	case Location::RIGHT: {
		const GEOLIB::Point & v (*(getPoint(k)));
		const GEOLIB::Point & w (*(getPoint(k+1)));
		if (w[1] < pnt[1] && pnt[1] <= v[1]) return EdgeType::CROSSING;
		else return EdgeType::INESSENTIAL;
		break;
	}
	case Location::BETWEEN:
	case Location::SOURCE:
	case Location::DESTINATION:
		return EdgeType::TOUCHING;
	default:
		return EdgeType::INESSENTIAL;
	}
}

void Polygon::calculateAxisAlignedBoundingBox ()
{
	size_t n_nodes (getNumberOfPoints());
	for (size_t k(0); k<n_nodes; k++) {
		_aabb.update ((*(getPoint(k))));
	}
}

void Polygon::ensureCWOrientation ()
{
	// *** pre processing: rotate points to xy-plan
	// *** copy points to vector - last point is identical to the first
	size_t n_pnts (this->getNumberOfPoints()-1);
	std::vector<GEOLIB::Point*> tmp_polygon_pnts;
	for (size_t k(0); k < n_pnts; k++) {
		tmp_polygon_pnts.push_back (new GEOLIB::Point (*(this->getPoint(k))));
	}

	// *** calculate supporting plane (plane normal and
	MathLib::Vector plane_normal;
	double d;
	MathLib::getNewellPlane(tmp_polygon_pnts, plane_normal, d);

	// *** rotate if necessary
	double tol (sqrt(std::numeric_limits<double>::min()));
	if (fabs(plane_normal[0]) > tol || fabs(plane_normal[1]) > tol) {
		// rotate copied points into x-y-plane
		MathLib::rotatePointsToXY(plane_normal, tmp_polygon_pnts);
	}

	for (size_t k(0); k<tmp_polygon_pnts.size(); k++) {
		(*(tmp_polygon_pnts[k]))[2] = 0.0; // should be -= d but there are numerical errors
	}

	// *** get the left most upper point
	size_t min_x_max_y_idx (0);	// for orientation check
	for (size_t k(0); k<n_pnts; k++) {
		if ((*(tmp_polygon_pnts[k]))[0] <= (*(tmp_polygon_pnts[min_x_max_y_idx]))[0]) {
			if ((*(tmp_polygon_pnts[k]))[0] < (*(tmp_polygon_pnts[min_x_max_y_idx]))[0]) {
				min_x_max_y_idx = k;
			} else {
				if ((*(tmp_polygon_pnts[k]))[1] > (*(tmp_polygon_pnts[min_x_max_y_idx]))[1]) {
					min_x_max_y_idx = k;
				}
			}
		}
	}
	// *** determine orientation
	MathLib::Orientation orient;
	if (0 < min_x_max_y_idx && min_x_max_y_idx < n_pnts-2) {
		orient = MathLib::getOrientation (
			tmp_polygon_pnts[min_x_max_y_idx-1],
			tmp_polygon_pnts[min_x_max_y_idx],
			tmp_polygon_pnts[min_x_max_y_idx+1]);
	} else {
		if (0 == min_x_max_y_idx) {
			orient = MathLib::getOrientation (
					tmp_polygon_pnts[n_pnts-1],
					tmp_polygon_pnts[0],
					tmp_polygon_pnts[1]);
		} else {
			orient = MathLib::getOrientation (
					tmp_polygon_pnts[n_pnts-2],
					tmp_polygon_pnts[n_pnts-1],
					tmp_polygon_pnts[0]);
		}
	}

	if (orient == MathLib::CCW) {
		// switch orientation
		size_t tmp_n_pnts (n_pnts);
		tmp_n_pnts++; // include last point of polygon (which is identical to the first)
		for (size_t k(0); k<tmp_n_pnts/2; k++) {
			BaseLib::swap (_ply_pnt_ids[k], _ply_pnt_ids[tmp_n_pnts-1-k]);
		}
	}

	for (size_t k(0); k<n_pnts; k++) {
		delete tmp_polygon_pnts[k];
	}
}

void Polygon::splitPolygonAtIntersection (std::list<Polygon*>::iterator polygon_it)
{
	size_t idx0 (0), idx1 (0);
	while (polygon_it != _simple_polygon_list.end()) {
		GEOLIB::Point *intersection_pnt (new GEOLIB::Point);
		bool is_simple (!MathLib::lineSegmentsIntersect (*polygon_it, idx0, idx1, *intersection_pnt));
		if (!is_simple) {
			// adding intersection point to pnt_vec
			size_t intersection_pnt_id (_ply_pnts.size());
			const_cast<std::vector<Point*>& >(_ply_pnts).push_back (intersection_pnt);

			// split Polygon
			if (idx0 > idx1) BaseLib::swap (idx0, idx1);

			GEOLIB::Polygon* polygon0 (new GEOLIB::Polygon((*polygon_it)->getPointsVec(), false));
			for (size_t k(0); k<=idx0; k++) polygon0->addPoint ((*polygon_it)->getPointID (k));
			polygon0->addPoint (intersection_pnt_id);
			for (size_t k(idx1+1); k<(*polygon_it)->getNumberOfPoints(); k++)
				polygon0->addPoint ((*polygon_it)->getPointID (k));
			if (! polygon0->initialise()) {
				std::cerr << "ERROR in Polygon::splitPolygonAtIntersection polygon0" << std::endl;
				exit (1);
			}

			GEOLIB::Polygon* polygon1 (new GEOLIB::Polygon((*polygon_it)->getPointsVec(), false));
			polygon1->addPoint (intersection_pnt_id);
			for (size_t k(idx0+1); k<=idx1; k++)
				polygon1->addPoint ((*polygon_it)->getPointID (k));
			polygon1->addPoint (intersection_pnt_id);
			if (! polygon1->initialise()) {
				std::cerr << "ERROR in Polygon::splitPolygonAtIntersection polygon1" << std::endl;
				exit (1);
			}

			// remove original polyline and add two new polylines
			std::list<GEOLIB::Polygon*>::iterator polygon0_it, polygon1_it;
			polygon_it = _simple_polygon_list.erase (polygon_it);
			polygon1_it = _simple_polygon_list.insert (polygon_it, polygon1);
			polygon0_it = _simple_polygon_list.insert (polygon1_it, polygon0);

			splitPolygonAtIntersection (polygon0_it);
			splitPolygonAtIntersection (polygon1_it);
		} else {
			delete intersection_pnt;
		}
		++polygon_it;
	}
}

void Polygon::splitPolygonAtPoint (std::list<GEOLIB::Polygon*>::iterator polygon_it)
{
	size_t n ((*polygon_it)->getNumberOfPoints()-1), idx0 (0), idx1(0);
	size_t *id_vec (new size_t[n]), *perm (new size_t[n]);
	for (size_t k(0); k<n; k++) {
		id_vec[k] = (*polygon_it)->getPointID (k);
		perm[k] = k;
	}

	BaseLib::quicksort (id_vec, 0, n, perm);

	for (size_t k(0); k<n-1; k++) {
		if (id_vec[k] == id_vec[k+1]) {
			idx0 = perm[k];
			idx1 = perm[k+1];
			delete [] perm;
			delete [] id_vec;

			if (idx0 > idx1) BaseLib::swap (idx0, idx1);

			// create two closed polylines
			GEOLIB::Polygon* polygon0 (new GEOLIB::Polygon((*polygon_it)->getPointsVec()));
			for (size_t k(0); k<=idx0; k++)
				polygon0->addPoint ((*polygon_it)->getPointID (k));
			for (size_t k(idx1+1); k<(*polygon_it)->getNumberOfPoints(); k++)
				polygon0->addPoint ((*polygon_it)->getPointID (k));
			polygon0->initialise();

			GEOLIB::Polygon* polygon1 (new GEOLIB::Polygon((*polygon_it)->getPointsVec()));
			for (size_t k(idx0); k<=idx1; k++)
				polygon1->addPoint ((*polygon_it)->getPointID (k));
			polygon1->initialise();

			// remove original polygon and add two new polygons
			std::list<GEOLIB::Polygon*>::iterator polygon0_it, polygon1_it;
			polygon1_it = _simple_polygon_list.insert (_simple_polygon_list.erase (polygon_it), polygon1);
			polygon0_it = _simple_polygon_list.insert (polygon1_it, polygon0);

			splitPolygonAtPoint (polygon0_it);
			splitPolygonAtPoint (polygon1_it);

			return;
		}
	}
	delete [] perm;
	delete [] id_vec;
}

GEOLIB::Polygon* createPolygonFromCircle (GEOLIB::Point const& middle_pnt, double radius,
		std::vector<GEOLIB::Point*> & pnts, size_t resolution)
{
	const size_t off_set (pnts.size());
	// create points
	double angle (2.0 * M_PI / resolution);
	for (size_t k(0); k<resolution; k++) {
		GEOLIB::Point *pnt (new GEOLIB::Point(middle_pnt.getData()));
		(*pnt)[0] += radius * cos (k*angle);
		(*pnt)[1] += radius * sin (k*angle);
		pnts.push_back (pnt);
	}

	// create polygon
	GEOLIB::Polygon* polygon (new GEOLIB::Polygon (pnts, false));
	for (size_t k(0); k<resolution; k++) {
		polygon->addPoint (k+off_set);
	}
	polygon->addPoint (off_set);

	return polygon;
}


} // end namespace GEOLIB
