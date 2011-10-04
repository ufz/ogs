/*
 * PointVec.cpp
 *
 *  Created on: Jun 11, 2010
 *      Author: TF
 */

#include "PointVec.h"
#include "BruteForceClosestPair.h"

namespace GEOLIB {

PointVec::PointVec (const std::string& name, std::vector<Point*>* points,
		std::map<std::string, size_t>* name_id_map, PointType type) :
			_pnt_vec (points), _name_id_map (name_id_map), _type(type), _name (name),
			_sqr_shortest_dist (std::numeric_limits<double>::max())
{
	assert (_pnt_vec);
	size_t number_of_all_input_pnts (_pnt_vec->size());

	makePntsUnique (_pnt_vec, _pnt_id_map);

	if (number_of_all_input_pnts - _pnt_vec->size() > 0)
		std::cerr << "WARNING: there are " << number_of_all_input_pnts - _pnt_vec->size() << " double points" << std::endl;
//	calculateShortestDistance ();
	calculateAxisAlignedBoundingBox ();
}

PointVec::~PointVec ()
{
	for (size_t k(0); k<_pnt_vec->size(); k++) {
		delete (*_pnt_vec)[k];
		(*_pnt_vec)[k] = NULL;
	}
	delete _pnt_vec;
	delete _name_id_map;
}

size_t PointVec::push_back (Point *pnt)
{
	_pnt_id_map.push_back (uniqueInsert(pnt));
	return _pnt_id_map[_pnt_id_map.size()-1];
}

void PointVec::push_back (Point *pnt, const std::string& name)
{
	std::map<std::string,size_t>::const_iterator it (_name_id_map->find (name));
	if (it != _name_id_map->end()) {
		std::cerr << "ERROR: PointVec::push_back (): two points with the same name" << std::endl;
		return;
	}

	size_t id (uniqueInsert (pnt));
	_pnt_id_map.push_back (id);
	(*_name_id_map)[name] = id;
}

size_t PointVec::uniqueInsert (Point* pnt)
{
	size_t n (_pnt_vec->size()), k;
	const double eps (std::numeric_limits<double>::epsilon());
	for (k=0; k<n; k++) {
		if (fabs((*((*_pnt_vec)[k]))[0]-(*pnt)[0]) < eps
				&&  fabs( (*((*_pnt_vec)[k]))[1]-(*pnt)[1]) < eps
				&&  fabs( (*((*_pnt_vec)[k]))[2]-(*pnt)[2]) < eps) {
			break;
		}
	}

	if(k==n) {
		_pnt_vec->push_back (pnt);
		// update shortest distance and bounding box
		for (size_t i(0); i<n; i++) {
			double sqr_dist (MathLib::sqrDist((*_pnt_vec)[i], (*_pnt_vec)[n]));
			if (sqr_dist < _sqr_shortest_dist)
				_sqr_shortest_dist = sqr_dist;
			aabb.update (*((*_pnt_vec)[n]));
		}
		return n;
	}

	delete pnt;
	pnt = NULL;
	return k;
}

std::vector<Point*> * PointVec::filterStations(const std::vector<PropertyBounds> &bounds) const
{
	std::vector<Point*> *tmpStations (new std::vector<Point*>);
	size_t size (_pnt_vec->size());
	for (size_t i=0; i<size; i++) {
		if (static_cast<Station*>((*_pnt_vec)[i])->inSelection(bounds)) tmpStations->push_back((*_pnt_vec)[i]);
	}
	return tmpStations;
}

bool PointVec::getElementIDByName (const std::string& name, size_t &id) const
{
	std::map<std::string,size_t>::const_iterator it (_name_id_map->find (name));

	if (it != _name_id_map->end()) {
		id = it->second;
		return true;
	} else return false;
}

const Point* PointVec::getElementByName (const std::string& name) const
{
	size_t id;
	bool ret (getElementIDByName (name, id));
	if (ret) {
		return (*_pnt_vec)[id];
	} else {
		return NULL;
	}
}

bool PointVec::getNameOfElement (const Point* data, std::string& name) const
{
	for (size_t k(0); k<_pnt_vec->size(); k++) {
		if ((*_pnt_vec)[k] == data) {
			return getNameOfElementByID (k, name);
		}
	}
	return false;
}

bool PointVec::getNameOfElementByID (size_t id, std::string& element_name) const
{
	if (! _name_id_map) return false;
	// search in map for id
	std::map<std::string,size_t>::const_iterator it (_name_id_map->begin());
	while (it != _name_id_map->end()) {
		if (it->second == id) {
			element_name = it->first;
			return true;
		}
		it++;
	}
	return false;
}

double PointVec::getShortestPointDistance () const
{
	return sqrt (_sqr_shortest_dist);
}

void PointVec::makePntsUnique (std::vector<GEOLIB::Point*>* pnt_vec, std::vector<size_t> &pnt_id_map)
{
	size_t n_pnts_in_file (pnt_vec->size());
	std::vector<size_t> perm;
	pnt_id_map.reserve (n_pnts_in_file);
	for (size_t k(0); k<n_pnts_in_file; k++) {
		perm.push_back (k);
		pnt_id_map.push_back(k);
	}

	// sort the points
	Quicksort<GEOLIB::Point*> (*pnt_vec, 0, n_pnts_in_file, perm);

	// unfortunately quicksort is not stable -
	// sort identical points by id - to make sorting stable
	double eps (sqrt(std::numeric_limits<double>::min()));
	// determine intervals with identical points to resort for stability of sorting
	std::vector<size_t> identical_pnts_interval;
	bool identical (false);
	for (size_t k=0; k<n_pnts_in_file-1; k++) {
		if ( fabs((*((*pnt_vec)[k+1]))[0]-(*((*pnt_vec)[k]))[0]) < eps
			&&  fabs( (*((*pnt_vec)[k+1]))[1]-(*((*pnt_vec)[k]))[1]) < eps
			&&  fabs( (*((*pnt_vec)[k+1]))[2]-(*((*pnt_vec)[k]))[2]) < eps) {
			// points are identical, sort by id
			if (!identical) identical_pnts_interval.push_back (k);
			identical = true;
		} else {
			if (identical) identical_pnts_interval.push_back (k+1);
			identical = false;
		}
	}
	if (identical) identical_pnts_interval.push_back (n_pnts_in_file);

	for (size_t i(0); i<identical_pnts_interval.size()/2; i++) {
		// bubble sort by id
		size_t beg (identical_pnts_interval[2*i]);
		size_t end (identical_pnts_interval[2*i+1]);
		for (size_t j (beg); j<end; j++) {
			for (size_t k (beg); k<end-1; k++) {
				if (perm[k] > perm[k+1]) std::swap (perm[k], perm[k+1]);
			}
		}
	}

	// check if there are identical points
	for (size_t k=0; k<n_pnts_in_file-1; k++) {
		if ( fabs((*((*pnt_vec)[k+1]))[0]-(*((*pnt_vec)[k]))[0]) < eps
				&&  fabs( (*((*pnt_vec)[k+1]))[1]-(*((*pnt_vec)[k]))[1]) < eps
				&&  fabs( (*((*pnt_vec)[k+1]))[2]-(*((*pnt_vec)[k]))[2]) < eps) {
			pnt_id_map[perm[k+1]] = pnt_id_map[perm[k]];
		}
	}

	// reverse permutation
	Quicksort<GEOLIB::Point*> (perm, 0, n_pnts_in_file, *pnt_vec);

	// remove the second, third, ... occurrence from vector
	for (size_t k(0); k<n_pnts_in_file; k++) {
		if (pnt_id_map[k] < k) {
			delete (*pnt_vec)[k];
			(*pnt_vec)[k] = NULL;
		}
	}
	// remove NULL-ptr from vector
	for (std::vector<GEOLIB::Point*>::iterator it(pnt_vec->begin()); it != pnt_vec->end(); ) {
		if (*it == NULL) {
			it = pnt_vec->erase (it);
		}
		else it++;
	}

	// renumber id-mapping
	size_t cnt (0);
	for (size_t k(0); k<n_pnts_in_file; k++) {
		if (pnt_id_map[k] == k) { // point not removed, if necessary: id change
			pnt_id_map[k] = cnt;
			cnt++;
		} else {
			pnt_id_map[k] = pnt_id_map[pnt_id_map[k]];
		}
	}

	// KR correct renumbering of indices
//	size_t cnt(0);
//	std::map<size_t, size_t> reg_ids;
//	for (size_t k(0); k < n_pnts_in_file; k++) {
//		if (pnt_id_map[k] == k) {
//			reg_ids.insert(std::pair<size_t, size_t>(k, cnt));
//			cnt++;
//		} else reg_ids.insert(std::pair<size_t, size_t>(k, reg_ids[pnt_id_map[k]]));
//	}
//	for (size_t k(0); k < n_pnts_in_file; k++)
//		pnt_id_map[k] = reg_ids[k];
}

void PointVec::calculateShortestDistance ()
{
	size_t i, j;
	BruteForceClosestPair (*_pnt_vec, i, j);
	_sqr_shortest_dist = MathLib::sqrDist ((*_pnt_vec)[i], (*_pnt_vec)[j]);
}

void PointVec::calculateAxisAlignedBoundingBox ()
{
	const size_t n_pnts (_pnt_vec->size());
	for (size_t i(0); i<n_pnts; i++) {
		aabb.update (*(*_pnt_vec)[i]);
	}
}

} // end namespace
