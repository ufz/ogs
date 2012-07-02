/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Polyline.cpp
 *
 * Created on 2010-06-21 by Thomas Fischer
 */

// Base
#include "swap.h"

#include "Polyline.h"

namespace GeoLib {

Polyline::Polyline(const std::vector<Point*>& pnt_vec) :
	GeoObject(), _ply_pnts(pnt_vec)
{
	_length.push_back (0.0);
}

Polyline::Polyline(const Polyline& ply) :
	GeoObject(), _ply_pnts (ply._ply_pnts)
{
	for (size_t k(0); k<ply.getNumberOfPoints(); ++k) {
		_ply_pnt_ids.push_back (ply.getPointID (k));
	}

	if (ply.getNumberOfPoints() > 0) {
		for (size_t k(0); k<ply.getNumberOfPoints(); ++k) {
			_length.push_back (ply.getLength (k));
		}
	}
}

void Polyline::write(std::ostream &os) const
{
	size_t size(_ply_pnt_ids.size());
	for (size_t k(0); k < size; k++) {
		os << *(_ply_pnts[_ply_pnt_ids[k]]) << std::endl;
	}
}

void Polyline::addPoint(size_t point_id)
{
	assert(point_id < _ply_pnts.size());
	size_t n_pnts (_ply_pnt_ids.size());
	_ply_pnt_ids.push_back(point_id);

	if (n_pnts > 0) {
		double act_dist (sqrt(MathLib::sqrDist (_ply_pnts[_ply_pnt_ids[n_pnts-1]], _ply_pnts[point_id])));
		double dist_until_now (0.0);
		if (n_pnts > 1)
			dist_until_now = _length[n_pnts - 1];

		_length.push_back (dist_until_now + act_dist);
	}
}

size_t Polyline::getNumberOfPoints() const
{
	return _ply_pnt_ids.size();
}

bool Polyline::isClosed() const
{
	if (_ply_pnt_ids.front() == _ply_pnt_ids.back())
		return true;
	else
		return false;
}

size_t Polyline::getPointID(size_t i) const
{
	assert(i < _ply_pnt_ids.size());
	return _ply_pnt_ids[i];
}

void Polyline::setPointID(size_t idx, size_t id)
{
	assert(idx < _ply_pnt_ids.size());
	_ply_pnt_ids[idx] = id;
}

const Point* Polyline::operator[](size_t i) const
{
	assert(i < _ply_pnt_ids.size());
	return _ply_pnts[_ply_pnt_ids[i]];
}

const Point* Polyline::getPoint(size_t i) const
{
	assert(i < _ply_pnt_ids.size());
 	return _ply_pnts[_ply_pnt_ids[i]];
}

std::vector<Point*> const& Polyline::getPointsVec () const
{
	return _ply_pnts;
}

double Polyline::getLength (size_t k) const
{
	assert(k < _length.size());
	return _length[k];
}

const std::vector<double>& Polyline::getLengthVec () const
{
	return _length;
}


Polyline* Polyline::constructPolylineFromSegments(const std::vector<Polyline*> &ply_vec, double prox)
{
	size_t nLines = ply_vec.size();

	Polyline* new_ply = new Polyline(*ply_vec[0]);
	std::vector<GeoLib::Point*> pnt_vec(new_ply->getPointsVec());

	std::vector<Polyline*> local_ply_vec;
	for (size_t i = 1; i < nLines; i++) {
		local_ply_vec.push_back(ply_vec[i]);
	}

	while (!local_ply_vec.empty()) {
		bool ply_found(false);
		prox *= prox; // square distance once to save time later
		for (std::vector<Polyline*>::iterator it=local_ply_vec.begin(); it!=local_ply_vec.end(); ++it)
		{
			if (pnt_vec == (*it)->getPointsVec())
			{
				size_t nPoints((*it)->getNumberOfPoints());

				//if (new_ply->getPointID(0) == (*it)->getPointID(0))
				if (pointsAreIdentical(pnt_vec, new_ply->getPointID(0), (*it)->getPointID(0), prox))
				{
					Polyline* tmp = new Polyline((*it)->getPointsVec());
					for (size_t k = 0; k < nPoints; k++)
						tmp->addPoint((*it)->getPointID(nPoints - k - 1));

					size_t new_ply_size(new_ply->getNumberOfPoints());
					for (size_t k = 1; k < new_ply_size; k++)
						tmp->addPoint(new_ply->getPointID(k));
					delete new_ply;
					new_ply = tmp;
					ply_found = true;
				}
				//else if (new_ply->getPointID(0) == (*it)->getPointID(nPoints-1))
				else if (pointsAreIdentical(pnt_vec, new_ply->getPointID(0), (*it)->getPointID(nPoints-1), prox))
				{
					Polyline* tmp = new Polyline(**it);
					size_t new_ply_size(new_ply->getNumberOfPoints());
					for (size_t k = 1; k < new_ply_size; k++)
						tmp->addPoint(new_ply->getPointID(k));
					delete new_ply;
					new_ply = tmp;
					ply_found = true;
				}
				//else if (new_ply->getPointID(new_ply->getNumberOfPoints()-1) == (*it)->getPointID(0))
				else if (pointsAreIdentical(pnt_vec, new_ply->getPointID(new_ply->getNumberOfPoints()-1), (*it)->getPointID(0), prox))
				{
					for (size_t k=1; k<nPoints; k++)
						new_ply->addPoint((*it)->getPointID(k));
					ply_found = true;
				}
				//else if (new_ply->getPointID(new_ply->getNumberOfPoints()-1) == (*it)->getPointID(nPoints-1))
				else if (pointsAreIdentical(pnt_vec, new_ply->getPointID(new_ply->getNumberOfPoints()-1), (*it)->getPointID(nPoints-1), prox))
				{
					for (size_t k=1; k<nPoints; k++)
						new_ply->addPoint((*it)->getPointID(nPoints-k-1));
					ply_found = true;
				}
				if (ply_found) {
					local_ply_vec.erase(it);
					break;
				}
			} else
				std::cout
						<< "Error in Polyline::contructPolylineFromSegments() - Line segments use different point vectors..."
						<< std::endl;
		}

		if (!ply_found) {
			std::cout
					<< "Error in Polyline::contructPolylineFromSegments() - Not all segments are connected..."
					<< std::endl;
			new_ply = NULL;
			break;
		}
	}
	return new_ply;
}

bool Polyline::pointsAreIdentical(const std::vector<Point*> &pnt_vec, size_t i, size_t j, double prox)
{
	if (i==j) return true;
	return (MathLib::checkDistance( *pnt_vec[i], *pnt_vec[j], prox ));
}

Polyline* Polyline::closePolyline(const Polyline& ply)
{
	if (ply.getNumberOfPoints()>2)
	{
		Polyline* new_ply = new Polyline(ply);
		if (ply.isClosed()) return new_ply;
		new_ply->addPoint(new_ply->getPointID(0));
		return new_ply;
	}
	std::cout << "Error in Polyline::closePolyline() - Input polyline needs to be composed of at least three points..." << std::endl;
	return NULL;
}

Location::type Polyline::getLocationOfPoint (size_t k, GeoLib::Point const & pnt) const
{
	assert (k<_ply_pnt_ids.size()-1);

	GeoLib::Point const& source (*(_ply_pnts[_ply_pnt_ids[k]]));
	GeoLib::Point const& dest (*(_ply_pnts[_ply_pnt_ids[k+1]]));
	GeoLib::Point a (dest[0]-source[0], dest[1]-source[1], dest[2]-source[2]); // vector
	GeoLib::Point b (pnt[0]-source[0], pnt[1]-source[1], pnt[2]-source[2]); // vector

	double det_2x2 (a[0]*b[1] - a[1]*b[0]);

	if (det_2x2 > std::numeric_limits<double>::epsilon()) return Location::LEFT;
	if (std::numeric_limits<double>::epsilon() < fabs(det_2x2)) return Location::RIGHT;
	if (a[0]*b[0] < 0.0 || a[1]*b[1] < 0.0) return Location::BEHIND;
	if (MathLib::sqrNrm2(&a) < MathLib::sqrNrm2(&b)) return Location::BEYOND;
	if (MathLib::sqrDist (&pnt, _ply_pnts[_ply_pnt_ids[k]]) < sqrt(std::numeric_limits<double>::min()))
		return Location::SOURCE;
	if (MathLib::sqrDist (&pnt, _ply_pnts[_ply_pnt_ids[k+1]]) < sqrt(std::numeric_limits<double>::min()))
		return Location::DESTINATION;
	return Location::BETWEEN;
}

std::ostream& operator<< (std::ostream &os, const Polyline &pl)
{
	pl.write (os);
	return os;
}

bool containsEdge (const Polyline& ply, size_t id0, size_t id1)
{
	if (id0 == id1) {
		std::cerr << "no valid edge id0 == id1 == " << id0 << std::endl;
		return false;
	}
	if (id0 > id1) BaseLib::swap (id0,id1);
	const size_t n (ply.getNumberOfPoints() - 1);
	for (size_t k(0); k<n; k++) {
		size_t ply_pnt0 (ply.getPointID (k));
		size_t ply_pnt1 (ply.getPointID (k+1));
		if (ply_pnt0 > ply_pnt1)
			BaseLib::swap (ply_pnt0, ply_pnt1);
		if (ply_pnt0 == id0 && ply_pnt1 == id1)
			return true;
	}
	return false;
}


} // end namespace GeoLib
