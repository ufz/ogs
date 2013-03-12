/**
 * \file
 * \author Thomas Fischeror
 * \date   2010-06-11
 * \brief  Implementation of the PointVec class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ThirdParty/logog
#include "logog/include/logog.hpp"

// GeoLib
#include "BruteForceClosestPair.h"
#include "PointVec.h"
#include "PointWithID.h"

// BaseLib
#include "quicksort.h"

// MathLib
#include "MathTools.h"

namespace GeoLib
{
PointVec::PointVec (const std::string& name, std::vector<Point*>* points,
                    std::map<std::string, std::size_t>* name_id_map, PointType type, double rel_eps) :
	TemplateVec<Point> (name, points, name_id_map),
	_type(type), _sqr_shortest_dist (std::numeric_limits<double>::max()),
	_aabb(points->begin(), points->end())
{
	assert (_data_vec);
	std::size_t number_of_all_input_pnts (_data_vec->size());

	rel_eps *= sqrt(MathLib::sqrDist (&(_aabb.getMinPoint()),&(_aabb.getMaxPoint())));
	makePntsUnique (_data_vec, _pnt_id_map, rel_eps);

	if (number_of_all_input_pnts - _data_vec->size() > 0)
		WARN("PointVec::PointVec(): there are %d double points.",
		     number_of_all_input_pnts - _data_vec->size());
}

PointVec::~PointVec ()
{}

std::size_t PointVec::push_back (Point* pnt)
{
	_pnt_id_map.push_back (uniqueInsert(pnt));
	return _pnt_id_map[_pnt_id_map.size() - 1];
}

void PointVec::push_back (Point* pnt, std::string const*const name)
{
	if (name == NULL) {
		_pnt_id_map.push_back (uniqueInsert(pnt));
		return;
	}

	std::map<std::string,std::size_t>::const_iterator it (_name_id_map->find (*name));
	if (it != _name_id_map->end()) {
		WARN("PointVec::push_back(): two points share the name %s.", name->c_str());
		return;
	}

	std::size_t id (uniqueInsert (pnt));
	_pnt_id_map.push_back (id);
	(*_name_id_map)[*name] = id;
}

std::size_t PointVec::uniqueInsert (Point* pnt)
{
	std::size_t n (_data_vec->size()), k;
	const double eps (std::numeric_limits<double>::epsilon());
	for (k = 0; k < n; k++)
		if (fabs((*((*_data_vec)[k]))[0] - (*pnt)[0]) < eps
		    &&  fabs( (*((*_data_vec)[k]))[1] - (*pnt)[1]) < eps
		    &&  fabs( (*((*_data_vec)[k]))[2] - (*pnt)[2]) < eps)
			break;

	if(k == n) {
		_data_vec->push_back (pnt);
		// update bounding box
		_aabb.update (*((*_data_vec)[n]));
		// update shortest distance
		for (std::size_t i(0); i < n; i++) {
			double sqr_dist (MathLib::sqrDist((*_data_vec)[i], (*_data_vec)[n]));
			if (sqr_dist < _sqr_shortest_dist)
				_sqr_shortest_dist = sqr_dist;
		}
		return n;
	}

	delete pnt;
	pnt = NULL;
	return k;
}

std::vector<Point*>* PointVec::filterStations(const std::vector<PropertyBounds> &bounds) const
{
	std::vector<Point*>* tmpStations (new std::vector<Point*>);
	std::size_t size (_data_vec->size());
	for (std::size_t i = 0; i < size; i++)
		if (static_cast<Station*>((*_data_vec)[i])->inSelection(bounds))
			tmpStations->push_back((*_data_vec)[i]);
	return tmpStations;
}

double PointVec::getShortestPointDistance () const
{
	return sqrt (_sqr_shortest_dist);
}

void PointVec::makePntsUnique (std::vector<GeoLib::Point*> *& pnt_vec,
                               std::vector<std::size_t> &pnt_id_map, double eps)
{
	std::size_t n_pnts_in_file(pnt_vec->size());
	std::vector<std::size_t> perm;
	pnt_id_map.reserve(n_pnts_in_file);
	for (std::size_t k(0); k < n_pnts_in_file; k++) {
		perm.push_back(k);
		pnt_id_map.push_back(k);
	}

	// sort the points
	BaseLib::Quicksort<GeoLib::Point*>(*pnt_vec, 0, n_pnts_in_file, perm);

	// unfortunately quicksort is not stable -
	// sort identical points by id - to make sorting stable
	// determine intervals with identical points to resort for stability of sorting
	std::vector<std::size_t> identical_pnts_interval;
	bool identical(false);
	for (std::size_t k = 0; k < n_pnts_in_file - 1; k++) {
		if (fabs((*((*pnt_vec)[k + 1]))[0] - (*((*pnt_vec)[k]))[0]) < eps
				&& fabs((*((*pnt_vec)[k + 1]))[1] - (*((*pnt_vec)[k]))[1]) < eps
				&& fabs((*((*pnt_vec)[k + 1]))[2] - (*((*pnt_vec)[k]))[2]) < eps) {
			// points are identical, sort by id
			if (!identical)
				identical_pnts_interval.push_back(k);
			identical = true;
		} else {
			if (identical)
				identical_pnts_interval.push_back(k + 1);
			identical = false;
		}
	}
	if (identical)
		identical_pnts_interval.push_back(n_pnts_in_file);

	for (std::size_t i(0); i < identical_pnts_interval.size() / 2; i++) {
		// bubble sort by id
		std::size_t beg(identical_pnts_interval[2 * i]);
		std::size_t end(identical_pnts_interval[2 * i + 1]);
		for (std::size_t j(beg); j < end; j++)
			for (std::size_t k(beg); k < end - 1; k++)
				if (perm[k] > perm[k + 1])
					std::swap(perm[k], perm[k + 1]);

	}

	// check if there are identical points
	for (std::size_t k = 0; k < n_pnts_in_file - 1; k++) {
		if (fabs((*((*pnt_vec)[k + 1]))[0] - (*((*pnt_vec)[k]))[0]) < eps
				&& fabs((*((*pnt_vec)[k + 1]))[1] - (*((*pnt_vec)[k]))[1]) < eps
				&& fabs((*((*pnt_vec)[k + 1]))[2] - (*((*pnt_vec)[k]))[2]) < eps) {
			pnt_id_map[perm[k + 1]] = pnt_id_map[perm[k]];
		}
	}

	// reverse permutation
	BaseLib::Quicksort<std::size_t, GeoLib::Point*>(perm, 0, n_pnts_in_file, *pnt_vec);

	// remove the second, third, ... occurrence from vector
	std::size_t cnt(0); // counts the points that are deleted
	for (std::size_t k(0); k < n_pnts_in_file; k++) {
		if (pnt_id_map[k] < k) {
			delete (*pnt_vec)[k];
			(*pnt_vec)[k] = nullptr;
			cnt++;
		}
	}

	std::vector<GeoLib::Point*>* tmp_pnt_vec(new std::vector<GeoLib::Point*>(n_pnts_in_file - cnt));
	std::remove_copy(pnt_vec->begin(), pnt_vec->end(), tmp_pnt_vec->begin(), nullptr);
	std::swap(tmp_pnt_vec, pnt_vec);
	delete tmp_pnt_vec;

	// renumber id-mapping
	cnt = 0;
	for (std::size_t k(0); k < n_pnts_in_file; k++) {
		if (pnt_id_map[k] == k) { // point not removed, if necessary: id change
			pnt_id_map[k] = cnt;
			cnt++;
		} else {
			pnt_id_map[k] = pnt_id_map[pnt_id_map[k]];
		}
	}
}

void PointVec::calculateShortestDistance ()
{
	std::size_t i, j;
	BruteForceClosestPair (*_data_vec, i, j);
	_sqr_shortest_dist = MathLib::sqrDist ((*_data_vec)[i], (*_data_vec)[j]);
}

std::vector<GeoLib::Point*>* PointVec::getSubset(const std::vector<std::size_t> &subset)
{
	std::vector<GeoLib::Point*> *new_points (new std::vector<GeoLib::Point*>(subset.size()));

	const std::size_t nPoints(subset.size());
	for (std::size_t i = 0; i < nPoints; i++)
		(*new_points)[i] = new GeoLib::PointWithID((*this->_data_vec)[subset[i]]->getCoords(), subset[i]);

	return new_points;
}
} // end namespace
