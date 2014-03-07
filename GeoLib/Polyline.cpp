/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-21
 * \brief  Implementation of the Polyline class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// STL
#include <algorithm>

// ThirdParty/logog
#include "logog/include/logog.hpp"

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
	for (std::size_t k(0); k < ply.getNumberOfPoints(); ++k)
		_ply_pnt_ids.push_back (ply.getPointID (k));

	if (ply.getNumberOfPoints() > 0)
		for (std::size_t k(0); k < ply.getNumberOfPoints(); ++k)
			_length.push_back (ply.getLength (k));
}

void Polyline::write(std::ostream &os) const
{
	std::size_t size(_ply_pnt_ids.size());
	for (std::size_t k(0); k < size; k++)
		os << *(_ply_pnts[_ply_pnt_ids[k]]) << "\n";
}

void Polyline::addPoint(std::size_t pnt_id)
{
	assert(pnt_id < _ply_pnts.size());
	std::size_t n_pnts (_ply_pnt_ids.size());

	// don't insert point if this would result in identical IDs for two adjacent points
	if (n_pnts > 0 && _ply_pnt_ids[n_pnts - 1] == pnt_id)
		return;

	_ply_pnt_ids.push_back(pnt_id);

	if (n_pnts > 0)
	{
		double act_dist (sqrt(MathLib::sqrDist (_ply_pnts[_ply_pnt_ids[n_pnts - 1]],
		                                        _ply_pnts[pnt_id])));
		double dist_until_now (0.0);
		if (n_pnts > 1)
			dist_until_now = _length[n_pnts - 1];

		_length.push_back (dist_until_now + act_dist);
	}
}

void Polyline::insertPoint(std::size_t pos, std::size_t pnt_id)
{
	assert(pnt_id < _ply_pnts.size());
	assert(pos <= _ply_pnt_ids.size());

	if (pos == _ply_pnt_ids.size()) {
		addPoint(pnt_id);
		return;
	}

	// check if inserting pnt_id would result in two identical IDs for adjacent points
	if (pos == 0 && pnt_id == _ply_pnt_ids[0]) {
		return;
	} else {
		if (pos != 0) {
			if (pos == (_ply_pnt_ids.size() - 1) && pnt_id == _ply_pnt_ids[pos]) {
				return;
			} else {
				if (pnt_id == _ply_pnt_ids[pos - 1] || pnt_id == _ply_pnt_ids[pos]) {
					return;
				}
			}
		}
	}

	std::vector<std::size_t>::iterator it(_ply_pnt_ids.begin() + pos);
	_ply_pnt_ids.insert(it, pnt_id);

	if (_ply_pnt_ids.size() > 1) {
		// update the _length vector
		if (pos == 0) {
			// insert at first position
			double act_dist(sqrt(MathLib::sqrDist(_ply_pnts[_ply_pnt_ids[1]],
			                                      _ply_pnts[pnt_id])));
			_length.insert(_length.begin() + 1, act_dist);
			const std::size_t s(_length.size());
			for (std::size_t k(2); k < s; k++)
				_length[k] += _length[1];
		} else {
			if (pos == _ply_pnt_ids.size() - 1) {
				// insert at last position
				double act_dist(sqrt(MathLib::sqrDist(
				                             _ply_pnts[_ply_pnt_ids[_ply_pnt_ids.size() - 2]],
				                             _ply_pnts[pnt_id])));
				double dist_until_now (0.0);
				if (_ply_pnt_ids.size() > 2)
					dist_until_now = _length[_ply_pnt_ids.size() - 2];

				_length.insert(_length.begin() + pos, dist_until_now + act_dist);
			} else {
				// insert at arbitrary position within the vector
				double dist_until_now (0.0);
				if (pos > 1)
					dist_until_now = _length[pos - 1];
				double len_seg0(sqrt(MathLib::sqrDist(
				                             _ply_pnts[_ply_pnt_ids[pos - 1]],
				                             _ply_pnts[pnt_id])));
				double len_seg1(sqrt(MathLib::sqrDist(
				                             _ply_pnts[_ply_pnt_ids[pos + 1]],
				                             _ply_pnts[pnt_id])));
				double update_dist(
				        len_seg0 + len_seg1 - (_length[pos] - dist_until_now));
				_length[pos] = dist_until_now + len_seg0;
				std::vector<double>::iterator it(_length.begin() + pos + 1);
				_length.insert(it, _length[pos] + len_seg1);
				for (it = _length.begin() + pos + 2; it != _length.end(); ++it)
					*it += update_dist;
			}
		}
	}
}

void Polyline::removePoint(std::size_t pos)
{
	if (pos >= _ply_pnt_ids.size())
		return;

	_ply_pnt_ids.erase(_ply_pnt_ids.begin() + pos);

	if (pos == _ply_pnt_ids.size())
	{
		_length.erase(_length.begin() + pos);
		return;
	}

	const std::size_t n_ply_pnt_ids(_ply_pnt_ids.size());
	if (pos == 0) {
		double seg_length(_length[0]);
		for (std::size_t k(0); k < n_ply_pnt_ids; k++)
			_length[k] = _length[k + 1] - seg_length;
		_length.pop_back();
	} else {
		const double len_seg0(_length[pos] - _length[pos - 1]);
		const double len_seg1(_length[pos + 1] - _length[pos]);
		_length.erase(_length.begin() + pos);
		const double len_new_seg(sqrt(MathLib::sqrDist(_ply_pnts[_ply_pnt_ids[pos - 1]],
		                                               _ply_pnts[_ply_pnt_ids[pos]])));
		double seg_length_diff(len_new_seg - len_seg0 - len_seg1);

		for (std::size_t k(pos); k < n_ply_pnt_ids; k++)
			_length[k] += seg_length_diff;
	}
}

std::size_t Polyline::getNumberOfPoints() const
{
	return _ply_pnt_ids.size();
}

bool Polyline::isClosed() const
{
	if (_ply_pnt_ids.size() < 3)
		return false;

	if (_ply_pnt_ids.front() == _ply_pnt_ids.back())
		return true;
	else
		return false;
}

bool Polyline::isPointIDInPolyline(std::size_t pnt_id) const
{
	return std::find(_ply_pnt_ids.begin(), _ply_pnt_ids.end(), pnt_id) != _ply_pnt_ids.end();
}

std::size_t Polyline::getPointID(std::size_t i) const
{
	assert(i < _ply_pnt_ids.size());
	return _ply_pnt_ids[i];
}

void Polyline::setPointID(std::size_t idx, std::size_t id)
{
	assert(idx < _ply_pnt_ids.size());
	_ply_pnt_ids[idx] = id;
}

const Point* Polyline::getPoint(std::size_t i) const
{
	assert(i < _ply_pnt_ids.size());
	return _ply_pnts[_ply_pnt_ids[i]];
}

std::vector<Point*> const& Polyline::getPointsVec () const
{
	return _ply_pnts;
}

double Polyline::getLength (std::size_t k) const
{
	assert(k < _length.size());
	return _length[k];
}

Polyline* Polyline::constructPolylineFromSegments(const std::vector<Polyline*> &ply_vec,
                                                  double prox)
{
	std::size_t nLines = ply_vec.size();

	Polyline* new_ply = new Polyline(*ply_vec[0]);
	std::vector<GeoLib::Point*> pnt_vec(new_ply->getPointsVec());

	std::vector<Polyline*> local_ply_vec;
	for (std::size_t i = 1; i < nLines; i++)
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
				std::size_t nPoints((*it)->getNumberOfPoints());

				//if (new_ply->getPointID(0) == (*it)->getPointID(0))
				if (pointsAreIdentical(pnt_vec, new_ply->getPointID(0),
				                       (*it)->getPointID(0), prox))
				{
					Polyline* tmp = new Polyline((*it)->getPointsVec());
					for (std::size_t k = 0; k < nPoints; k++)
						tmp->addPoint((*it)->getPointID(nPoints - k - 1));

					std::size_t new_ply_size(new_ply->getNumberOfPoints());
					for (std::size_t k = 1; k < new_ply_size; k++)
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
					std::size_t new_ply_size(new_ply->getNumberOfPoints());
					for (std::size_t k = 1; k < new_ply_size; k++)
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
					for (std::size_t k = 1; k < nPoints; k++)
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
					for (std::size_t k = 1; k < nPoints; k++)
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
				ERR("Error in Polyline::contructPolylineFromSegments() - Line segments use different point vectors.");
		}

		if (!ply_found)
		{
			ERR("Error in Polyline::contructPolylineFromSegments() - Not all segments are connected.");
			new_ply = NULL;
			break;
		}
	}
	return new_ply;
}

void Polyline::closePolyline()
{
	if (getNumberOfPoints() < 2) {
		ERR("Polyline::closePolyline(): Input polyline needs to be composed of at least three points.");
	}
	if (! isClosed()) {
		addPoint(getPointID(0));
	}
}

Location Polyline::getLocationOfPoint (std::size_t k, GeoLib::Point const & pnt) const
{
	assert (k < _ply_pnt_ids.size() - 1);

	GeoLib::Point const& source (*(_ply_pnts[_ply_pnt_ids[k]]));
	GeoLib::Point const& dest (*(_ply_pnts[_ply_pnt_ids[k + 1]]));
	long double a[2] = {dest[0] - source[0], dest[1] - source[1]}; // vector
	long double b[2] = {pnt[0] - source[0], pnt[1] - source[1]}; // vector

	long double det_2x2 (a[0] * b[1] - a[1] * b[0]);

	if (det_2x2 > std::numeric_limits<double>::epsilon())
		return Location::LEFT;
	if (std::numeric_limits<double>::epsilon() < std::abs(det_2x2))
		return Location::RIGHT;
	if (a[0] * b[0] < 0.0 || a[1] * b[1] < 0.0)
		return Location::BEHIND;
	if (a[0] * a[0] + a[1] * a[1] < b[0] * b[0] + b[1] * b[1])
		return Location::BEYOND;
	if (MathLib::sqrDist (&pnt,
	                      _ply_pnts[_ply_pnt_ids[k]]) < pow(std::numeric_limits<double>::epsilon(),2))
		return Location::SOURCE;
	if (MathLib::sqrDist (&pnt,
	                      _ply_pnts[_ply_pnt_ids[k + 1]]) <
	    sqrt(std::numeric_limits<double>::epsilon()))
		return Location::DESTINATION;
	return Location::BETWEEN;
}

void Polyline::updatePointIDs(const std::vector<std::size_t> &pnt_ids)
{
	for (auto it = this->_ply_pnt_ids.begin(); it!=this->_ply_pnt_ids.end();)
	{
		if (pnt_ids[*it] != *it)
		{
			if (it!=this->_ply_pnt_ids.begin() && (pnt_ids[*it] == pnt_ids[*(it-1)]))
				it = this->_ply_pnt_ids.erase(it);
			else
			{
				*it = pnt_ids[*it];
				++it;
			}
		}
		else
			++it;
	}
}

std::ostream& operator<< (std::ostream &os, const Polyline &pl)
{
	pl.write (os);
	return os;
}

bool containsEdge (const Polyline& ply, std::size_t id0, std::size_t id1)
{
	if (id0 == id1)
	{
		ERR("no valid edge id0 == id1 == %d.", id0);
		return false;
	}
	if (id0 > id1)
		std::swap (id0,id1);
	const std::size_t n (ply.getNumberOfPoints() - 1);
	for (std::size_t k(0); k < n; k++)
	{
		std::size_t ply_pnt0 (ply.getPointID (k));
		std::size_t ply_pnt1 (ply.getPointID (k + 1));
		if (ply_pnt0 > ply_pnt1)
			std::swap (ply_pnt0, ply_pnt1);
		if (ply_pnt0 == id0 && ply_pnt1 == id1)
			return true;
	}
	return false;
}

bool isLineSegmentIntersecting (const Polyline& ply,
                                GeoLib::Point const& s0,
                                GeoLib::Point const& s1)
{
	const std::size_t n (ply.getNumberOfPoints() - 1);
	bool intersect(false);
	GeoLib::Point intersection_pnt;
	for (std::size_t k(0); k < n && !intersect; k++)
		intersect = GeoLib::lineSegmentIntersect (*(ply.getPoint(k)), *(ply.getPoint(
		                                                                         k + 1)),
		                                           s0, s1, intersection_pnt);
	return intersect;
}

bool operator==(Polyline const& lhs, Polyline const& rhs)
{
	if (lhs.getNumberOfPoints() != rhs.getNumberOfPoints())
		return false;

	const std::size_t n(lhs.getNumberOfPoints());
	for (std::size_t k(0); k < n; k++)
		if (lhs.getPointID(k) != rhs.getPointID(k))
			return false;

	return true;
}

bool pointsAreIdentical(const std::vector<Point*> &pnt_vec,
                        std::size_t i,
                        std::size_t j,
                        double prox)
{
	if (i == j)
		return true;
	return MathLib::sqrDist(pnt_vec[i], pnt_vec[j]) < prox;
}
} // end namespace GeoLib
