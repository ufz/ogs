/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file Polyline.cpp
 *
 *  Created on 2010-06-21 by Thomas Fischer
 */

// Base

// GeoLib
#include "Polyline.h"

// MathLib
#include "AnalyticalGeometry.h"

namespace GeoLib
{
Polyline::Polyline(const std::vector<Point*>& pnt_vec) :
	GeoObject(), _ply_pnts(pnt_vec)
{
	_length.push_back (0.0);
}

Polyline::Polyline(const Polyline& ply) :
	GeoObject(), _ply_pnts (ply._ply_pnts)
{
	for (size_t k(0); k < ply.getNumberOfPoints(); ++k)
		_ply_pnt_ids.push_back (ply.getPointID (k));

	if (ply.getNumberOfPoints() > 0)
		for (size_t k(0); k < ply.getNumberOfPoints(); ++k)
			_length.push_back (ply.getLength (k));

}

void Polyline::write(std::ostream &os) const
{
	size_t size(_ply_pnt_ids.size());
	for (size_t k(0); k < size; k++)
		os << *(_ply_pnts[_ply_pnt_ids[k]]) << std::endl;
}

void Polyline::addPoint(size_t pnt_id)
{
	assert(pnt_id < _ply_pnts.size());
	size_t n_pnts (_ply_pnt_ids.size());

	// don't insert point if ID if this would result in identical IDs for two adjacent points
	if (n_pnts>0 && _ply_pnt_ids[n_pnts-1] == pnt_id) return;

	_ply_pnt_ids.push_back(pnt_id);

	if (n_pnts > 0) {
		double act_dist (sqrt(MathLib::sqrDist (_ply_pnts[_ply_pnt_ids[n_pnts - 1]],
		                                        _ply_pnts[pnt_id])));
		double dist_until_now (0.0);
		if (n_pnts > 1)
			dist_until_now = _length[n_pnts - 1];

		_length.push_back (dist_until_now + act_dist);
	}
}

void Polyline::insertPoint(size_t pos, size_t pnt_id)
{
	assert(pnt_id < _ply_pnts.size());
	assert(pos < _ply_pnt_ids.size());

	// check if inserting pnt_id would result in two identical IDs for adjacent points
	if (pos == 0 && pnt_id == _ply_pnt_ids[0]) return;
	else if (pos == (_ply_pnt_ids.size()-1) && pnt_id == _ply_pnt_ids[pos]) return;
	else if (pnt_id == _ply_pnt_ids[pos-1] || pnt_id == _ply_pnt_ids[pos]) return;

	std::vector<size_t>::iterator it(_ply_pnt_ids.begin() + pos);
	_ply_pnt_ids.insert(it, pnt_id);

	if (_ply_pnt_ids.size() > 1) {
		// update the _length vector
		if (pos == 0) {
			// insert at first position
			double act_dist(sqrt(MathLib::sqrDist(_ply_pnts[_ply_pnt_ids[1]], _ply_pnts[pnt_id])));
			_length.insert(_length.begin(), act_dist);
			const size_t s(_length.size());
			for (size_t k(1); k<s; k++) {
				_length[k] += _length[0];
			}
		} else {
			if (pos == _ply_pnt_ids.size()-1) {
				// insert at last position
				double act_dist(sqrt(MathLib::sqrDist(_ply_pnts[_ply_pnt_ids[_ply_pnt_ids.size()-2]], _ply_pnts[pnt_id])));
				double dist_until_now (0.0);
				if (_ply_pnt_ids.size() > 2)
					dist_until_now = _length[_ply_pnt_ids.size() - 2];

				_length.insert(_length.begin()+pos, dist_until_now + act_dist);
			} else {
				// insert at arbitrary position within the vector
				double dist_until_now (0.0);
				if (pos > 1) {
					dist_until_now = _length[pos-2];
				}
				double len_seg0(sqrt(MathLib::sqrDist(_ply_pnts[_ply_pnt_ids[pos-1]], _ply_pnts[pnt_id])));
				double len_seg1(sqrt(MathLib::sqrDist(_ply_pnts[_ply_pnt_ids[pos+1]], _ply_pnts[pnt_id])));
				_length[pos-1] = dist_until_now + len_seg0;
				std::vector<double>::iterator it(_length.begin()+pos);
				_length.insert(it, _length[pos-1]+len_seg1);
			}
		}
	}
}

void Polyline::removePoint(std::size_t pos)
{
	if (pos >= _ply_pnt_ids.size())
		return;

	_ply_pnt_ids.erase(_ply_pnt_ids.begin()+pos);

	if (pos == _ply_pnt_ids.size()-1) {
		_length.erase(_length.begin()+pos);
		return;
	}

	const size_t n_ply_pnt_ids(_ply_pnt_ids.size());
	if (pos == 0) {
		double seg_length(_length[0]);
		for (unsigned k(0); k<n_ply_pnt_ids; k++) {
			_length[k] = _length[k+1] - seg_length;
		}
		_length.pop_back();
	} else {
		const double len_seg0(_length[pos] - _length[pos-1]);
		const double len_seg1(_length[pos+1] - _length[pos]);
		_length.erase(_length.begin()+pos);
		const double len_new_seg(sqrt(MathLib::sqrDist(_ply_pnts[_ply_pnt_ids[pos-1]], _ply_pnts[pos])));
		double seg_length_diff(len_new_seg - len_seg0 - len_seg1);

		for (unsigned k(pos); k<n_ply_pnt_ids; k++) {
			_length[k] += seg_length_diff;
		}
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

bool Polyline::isPointIDInPolyline(size_t pnt_id) const
{
	const size_t n_ply_pnt_ids(_ply_pnt_ids.size());
	size_t k(0);
	while (k<n_ply_pnt_ids && _ply_pnt_ids[k] != pnt_id) {
		k++;
	}

	if (k == n_ply_pnt_ids) {
		return false;
	}
	return true;
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

Polyline* Polyline::constructPolylineFromSegments(const std::vector<Polyline*> &ply_vec,
                                                  double prox)
{
	size_t nLines = ply_vec.size();

	Polyline* new_ply = new Polyline(*ply_vec[0]);
	std::vector<GeoLib::Point*> pnt_vec(new_ply->getPointsVec());

	std::vector<Polyline*> local_ply_vec;
	for (size_t i = 1; i < nLines; i++)
		local_ply_vec.push_back(ply_vec[i]);

	while (!local_ply_vec.empty())
	{
		bool ply_found(false);
		prox *= prox; // square distance once to save time later
		for (std::vector<Polyline*>::iterator it = local_ply_vec.begin();
		     it != local_ply_vec.end(); ++it)
		{
			if (pnt_vec == (*it)->getPointsVec())
			{
				size_t nPoints((*it)->getNumberOfPoints());

				//if (new_ply->getPointID(0) == (*it)->getPointID(0))
				if (pointsAreIdentical(pnt_vec, new_ply->getPointID(0),
				                       (*it)->getPointID(0), prox))
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
				else if (pointsAreIdentical(pnt_vec, new_ply->getPointID(0),
				                            (*it)->getPointID(nPoints - 1), prox))
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
				else if (pointsAreIdentical(pnt_vec,
				                            new_ply->getPointID(new_ply->
				                                                getNumberOfPoints()
				                                                - 1),
				                            (*it)->getPointID(0), prox))
				{
					for (size_t k = 1; k < nPoints; k++)
						new_ply->addPoint((*it)->getPointID(k));
					ply_found = true;
				}
				//else if (new_ply->getPointID(new_ply->getNumberOfPoints()-1) == (*it)->getPointID(nPoints-1))
				else if (pointsAreIdentical(pnt_vec,
				                            new_ply->getPointID(new_ply->
				                                                getNumberOfPoints()
				                                                - 1),
				                            (*it)->getPointID(nPoints - 1), prox))
				{
					for (size_t k = 1; k < nPoints; k++)
						new_ply->addPoint((*it)->getPointID(nPoints - k - 1));
					ply_found = true;
				}
				if (ply_found)
				{
					local_ply_vec.erase(it);
					break;
				}
			}
			else
				std::cout
				<<
				"Error in Polyline::contructPolylineFromSegments() - Line segments use different point vectors..."
				<< std::endl;
		}

		if (!ply_found)
		{
			std::cout
			<<
			"Error in Polyline::contructPolylineFromSegments() - Not all segments are connected..."
			<< std::endl;
			new_ply = NULL;
			break;
		}
	}
	return new_ply;
}

bool Polyline::pointsAreIdentical(const std::vector<Point*> &pnt_vec,
                                  size_t i,
                                  size_t j,
                                  double prox)
{
	if (i == j)
		return true;
	return MathLib::checkDistance( *pnt_vec[i], *pnt_vec[j], prox );
}

Polyline* Polyline::closePolyline(const Polyline& ply)
{
	if (ply.getNumberOfPoints() > 2)
	{
		Polyline* new_ply = new Polyline(ply);
		if (ply.isClosed())
			return new_ply;
		new_ply->addPoint(new_ply->getPointID(0));
		return new_ply;
	}
	std::cout <<
	"Error in Polyline::closePolyline() - Input polyline needs to be composed of at least three points..."
	          << std::endl;
	return NULL;
}

Location::type Polyline::getLocationOfPoint (size_t k, GeoLib::Point const & pnt) const
{
	assert (k < _ply_pnt_ids.size() - 1);

	GeoLib::Point const& source (*(_ply_pnts[_ply_pnt_ids[k]]));
	GeoLib::Point const& dest (*(_ply_pnts[_ply_pnt_ids[k + 1]]));
	long double a[2] = {dest[0] - source[0], dest[1] - source[1]}; // vector
	long double b[2] = {pnt[0] - source[0], pnt[1] - source[1]}; // vector

	long double det_2x2 (a[0] * b[1] - a[1] * b[0]);

	if (det_2x2 > std::numeric_limits<double>::epsilon())
		return Location::LEFT;
	if (std::numeric_limits<double>::epsilon() < fabs(det_2x2))
		return Location::RIGHT;
	if (a[0] * b[0] < 0.0 || a[1] * b[1] < 0.0)
		return Location::BEHIND;
	if (a[0]*a[0]+a[1]*a[1] < b[0]*b[0]+b[1]*b[1])
		return Location::BEYOND;
	if (MathLib::sqrDist (&pnt, _ply_pnts[_ply_pnt_ids[k]]) < sqrt(std::numeric_limits<double>::min()))
		return Location::SOURCE;
	if (MathLib::sqrDist (&pnt, _ply_pnts[_ply_pnt_ids[k + 1]]) < sqrt(std::numeric_limits<double>::min()))
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
	if (id0 == id1)
	{
		std::cerr << "no valid edge id0 == id1 == " << id0 << std::endl;
		return false;
	}
	if (id0 > id1)
		std::swap (id0,id1);
	const size_t n (ply.getNumberOfPoints() - 1);
	for (size_t k(0); k < n; k++)
	{
		size_t ply_pnt0 (ply.getPointID (k));
		size_t ply_pnt1 (ply.getPointID (k + 1));
		if (ply_pnt0 > ply_pnt1)
			std::swap (ply_pnt0, ply_pnt1);
		if (ply_pnt0 == id0 && ply_pnt1 == id1)
			return true;
	}
	return false;
}

bool isLineSegmentIntersecting (const Polyline& ply, GeoLib::Point const& s0, GeoLib::Point const& s1)
{
	const size_t n (ply.getNumberOfPoints() - 1);
	bool intersect(false);
	GeoLib::Point intersection_pnt;
	for (size_t k(0); k < n && !intersect; k++) {
		intersect = MathLib::lineSegmentIntersect (*(ply.getPoint(k)), *(ply.getPoint(k+1)), s0, s1, intersection_pnt);
	}
	return intersect;
}

bool operator==(Polyline const& lhs, Polyline const& rhs)
{
	if (lhs.getNumberOfPoints() != rhs.getNumberOfPoints())
		return false;

	const size_t n(lhs.getNumberOfPoints());
	for (size_t k(0); k<n; k++) {
		if (lhs.getPointID(k) != rhs.getPointID(k))
			return false;
	}

	return true;
}

} // end namespace GeoLib
