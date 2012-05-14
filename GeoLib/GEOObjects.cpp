/*
 * GEOObjects.cpp
 *
 *  Created on: Jan 21, 2010
 *      Author: TF / KR
 */

#include "GEOObjects.h"
#include "StringTools.h"

#include <fstream>

namespace GEOLIB {

GEOObjects::GEOObjects()
{
}

GEOObjects::~GEOObjects()
{
	// delete surfaces
	for (size_t k(0); k < _sfc_vecs.size(); k++) {
		delete _sfc_vecs[k];
	}
	// delete polylines
	for (size_t k(0); k < _ply_vecs.size(); k++) {
		delete _ply_vecs[k];
	}
	// delete points
	for (size_t k(0); k < _pnt_vecs.size(); k++) {
		delete _pnt_vecs[k];
	}
}

void GEOObjects::addPointVec(std::vector<Point*> *points, std::string &name, std::map<std::string, size_t>* pnt_id_name_map)
{
	isUniquePointVecName(name);
	_pnt_vecs.push_back(new PointVec(name, points, pnt_id_name_map));
//	std::cout << "minimal dist between points: " << (_pnt_vecs[_pnt_vecs.size()-1])->getShortestPointDistance () << std::endl;
}

bool GEOObjects::appendPointVec(std::vector<Point*> const& new_points,
		std::string const &name, std::vector<size_t>* ids)
{
	// search vector
	size_t idx (0);
	bool nfound (true);
	for (idx=0; idx<_pnt_vecs.size() && nfound; idx++) {
		if ( (_pnt_vecs[idx]->getName()).compare (name) == 0 )
			nfound = false;
	}

	if (! nfound) {
		idx--;
		size_t n_new_pnts (new_points.size());
		// append points
		if (ids) {
			for (size_t k(0); k<n_new_pnts; k++) {
				ids->push_back (_pnt_vecs[idx]->push_back (new_points[k]));
			}
		} else {
			for (size_t k(0); k<n_new_pnts; k++) {
				_pnt_vecs[idx]->push_back (new_points[k]);
			}
		}
		return true;
	} else return false;
}

const std::vector<Point*>* GEOObjects::getPointVec(const std::string &name) const
{
	size_t size (_pnt_vecs.size());
	for (size_t i=0; i<size; i++)
	{
		if (_pnt_vecs[i]->getName().compare(name)==0)
			return _pnt_vecs[i]->getVector();
	}
	std::cout << "GEOObjects::getPointVec() - No entry found with name \"" << name << "\"." << std::endl;
	return NULL;
}

const PointVec* GEOObjects::getPointVecObj(const std::string &name) const
{
	size_t size (_pnt_vecs.size());
	for (size_t i=0; i<size; i++) {
		if (_pnt_vecs[i]->getName().compare(name)==0)
			return _pnt_vecs[i];
	}
	std::cout << "GEOObjects::getPointVec() - No entry found with name \"" << name << "\"." << std::endl;
	return NULL;
}

bool GEOObjects::removePointVec(const std::string &name)
{
	if (isPntVecUsed (name)) {
		std::cout << "GEOObjects::removePointVec() - There are still Polylines or Surfaces depending on these points." << std::endl;
		return false;
	}

	for (std::vector<PointVec*>::iterator it(_pnt_vecs.begin());
		it != _pnt_vecs.end(); it++) {
		if ((*it)->getName().compare(name)==0) {
			delete *it;
			_pnt_vecs.erase(it);
			return true;
		}
	}
	std::cout << "GEOObjects::removePointVec() - No entry found with name \"" << name << "." << std::endl;
	return false;
}


void GEOObjects::addStationVec(std::vector<Point*> *stations, std::string &name, const Color* const color)
{
	size_t size = stations->size();
	for (size_t i=0; i<size; i++) static_cast<Station*>((*stations)[i])->setColor(color);
	isUniquePointVecName(name);
	_pnt_vecs.push_back(new PointVec(name, stations, NULL, PointVec::STATION));
}


std::vector<Point*> *GEOObjects::filterStationVec(const std::string &name,
		const std::vector<PropertyBounds> &bounds)
{
	for (std::vector<PointVec*>::iterator it(_pnt_vecs.begin());
			it != _pnt_vecs.end(); it++) {
		if ((*it)->getName().compare(name) == 0 && (*it)->getType()
				== PointVec::STATION) {
			return (*it)->filterStations(bounds);
		}
	}
	std::cout << "GEOObjects:: filterStations() - No entry found with name \""
			<< name << "." << std::endl;
	return NULL;
}

const std::vector<Point*> *GEOObjects::getStationVec(const std::string &name) const
{
	for (std::vector<PointVec*>::const_iterator it(_pnt_vecs.begin());
		it != _pnt_vecs.end(); it++) {
		if ((*it)->getName().compare(name) == 0 && (*it)->getType()
				== PointVec::STATION)
			return (*it)->getVector();
	}
	std::cout << "GEOObjects::getStationVec() - No entry found with name \""
			<< name << "." << std::endl;
	return NULL;
}

void GEOObjects::addPolylineVec(std::vector<Polyline*> *lines,
		const std::string &name, std::map<std::string, size_t>* ply_names)
{
	for (std::vector<Polyline*>::iterator it (lines->begin());
		it != lines->end(); ) {
		if ((*it)->getNumberOfPoints() < 2) {
			std::vector<Polyline*>::iterator it_erase (it);
			it = lines->erase (it_erase);
		} else it++;
	}

	if (lines->empty()) return;

	_ply_vecs.push_back(new PolylineVec(name, lines, ply_names));
}

bool GEOObjects::appendPolylineVec(const std::vector<Polyline*> &polylines, const std::string &name)
{
	// search vector
	size_t idx (0);
	bool nfound (true);
	for (idx=0; idx<_ply_vecs.size() && nfound; idx++) {
		if ( (_ply_vecs[idx]->getName()).compare (name) == 0 ) nfound = false;
	}

	if (! nfound) {
		idx--;
		size_t n_plys (polylines.size());
		// append lines
		for (size_t k(0); k<n_plys; k++)
			_ply_vecs[idx]->push_back (polylines[k]);
		return true;
	} else return false;
}

const std::vector<Polyline*> *GEOObjects::getPolylineVec(const std::string &name) const
{
	size_t size (_ply_vecs.size());
	for (size_t i=0; i<size; i++)
	{
		if (_ply_vecs[i]->getName().compare(name)==0)
			return _ply_vecs[i]->getVector();
	}
#ifndef NDEBUG
	std::cout << "DEB: GEOObjects::getPolylineVec() - No entry found with name \"" << name << "." << std::endl;
#endif
	return NULL;
}

const PolylineVec* GEOObjects::getPolylineVecObj(const std::string &name) const
{
	size_t size (_ply_vecs.size());
	for (size_t i=0; i<size; i++) {
		if (_ply_vecs[i]->getName().compare(name)==0)
			return _ply_vecs[i];
	}
#ifndef NDEBUG
	std::cout << "DEB: GEOObjects::getPolylineVec() - No entry found with name \"" << name << "\"." << std::endl;
#endif
	return NULL;
}

bool GEOObjects::removePolylineVec(const std::string &name)
{
	for (std::vector<PolylineVec*>::iterator it = _ply_vecs.begin();
		it != _ply_vecs.end(); ++it)
	{
		if ((*it)->getName().compare(name) == 0) {
			delete *it;
			_ply_vecs.erase(it);
			return true;
		}
	}
#ifndef NDEBUG
	std::cout << "GEOObjects::removePolylineVec() - No entry found with name \""
			<< name << "\"." << std::endl;
#endif
	return false;
}

void GEOObjects::addSurfaceVec(std::vector<Surface*> *sfc, const std::string &name,
		std::map<std::string, size_t>* sfc_names)
{
	_sfc_vecs.push_back(new SurfaceVec(name, sfc, sfc_names));
}

bool GEOObjects::appendSurfaceVec(const std::vector<Surface*> &surfaces, const std::string &name)
{
	// search vector
	size_t idx (0);
	bool nfound (true);
	for (idx=0; idx<_sfc_vecs.size() && nfound; idx++) {
		if ( (_sfc_vecs[idx]->getName()).compare (name) == 0 ) nfound = false;
	}

	if (! nfound) {
		idx--;
		size_t n_sfcs (surfaces.size());
		// append surfaces
		for (size_t k(0); k<n_sfcs; k++)
			_sfc_vecs[idx]->push_back (surfaces[k]);
		return true;
	} else return false;
}

const std::vector<Surface*>* GEOObjects::getSurfaceVec(const std::string &name) const
{
	size_t size (_sfc_vecs.size());
	for (size_t i=0; i<size; i++)
	{
		if (_sfc_vecs[i]->getName().compare(name)==0)
			return _sfc_vecs[i]->getVector();
	}
	std::cout << "GEOObjects::getSurfaceVec() - No entry found with name \"" << name << "." << std::endl;
	return NULL;
}

bool GEOObjects::removeSurfaceVec(const std::string &name)
{
	for (std::vector<SurfaceVec*>::iterator it (_sfc_vecs.begin());
		it != _sfc_vecs.end(); it++) {
		if ((*it)->getName().compare (name) == 0) {
			delete *it;
			_sfc_vecs.erase (it);
			return true;
		}
	}
#ifndef NDEBUG
	std::cout << "GEOObjects::removeSurfaceVec() - No entry found with name \""
			<< name << "\"." << std::endl;
#endif
	return false;
}

const SurfaceVec* GEOObjects::getSurfaceVecObj(const std::string &name) const
{
	size_t size (_sfc_vecs.size());
	for (size_t i=0; i<size; i++) {
		if (_sfc_vecs[i]->getName().compare(name)==0)
			return _sfc_vecs[i];
	}
	std::cout << "GEOObjects::getSurfaceVec() - No entry found with name \"" << name << "\"." << std::endl;
	return NULL;
}

bool GEOObjects::isUniquePointVecName(std::string &name)
{
	int count=0;
	bool isUnique = false;
	std::string cpName;

	while (!isUnique)
	{
		isUnique = true;
		cpName = name;

		count++;
		// If the original name already exists we start to add numbers to name for
		// as long as it takes to make the name unique.
		if (count>1) cpName = cpName + "-" + number2str(count);

		for (size_t i=0; i<_pnt_vecs.size(); i++)
		{
			if ( cpName.compare(_pnt_vecs[i]->getName()) == 0 ) isUnique = false;
		}
	}

	// At this point cpName is a unique name and isUnique is true.
	// If cpName is not the original name, "name" is changed and isUnique is set to false,
	// indicating that a vector with the original name already exists.
	if (count>1)
	{
		isUnique = false;
		name = cpName;
	}
	return isUnique;
}

bool GEOObjects::isPntVecUsed (const std::string &name) const
{
	// search dependent data structures (Polyline)
	for (std::vector<PolylineVec*>::const_iterator it ( _ply_vecs.begin());	it != _ply_vecs.end(); it++)
	{
		std::string a = (*it)->getName();
		if (((*it)->getName()).compare(name) == 0)
			return true;
	}
	for (std::vector<SurfaceVec*>::const_iterator it ( _sfc_vecs.begin());	it != _sfc_vecs.end(); it++)
	{
		std::string a = (*it)->getName();
		if (((*it)->getName()).compare(name) == 0)
			return true;
	}

	return false;

}

void GEOObjects::getStationNames(std::vector<std::string>& names) const
{
	for (std::vector<PointVec*>::const_iterator it(_pnt_vecs.begin());	it != _pnt_vecs.end(); it++) {
		if ((*it)->getType() == PointVec::STATION)
			names.push_back((*it)->getName());
	}
}

void GEOObjects::getGeometryNames (std::vector<std::string>& names) const
{
	names.clear ();
	for (std::vector<PointVec*>::const_iterator it(_pnt_vecs.begin());	it != _pnt_vecs.end(); it++) {
		if ((*it)->getType() == PointVec::POINT)
			names.push_back((*it)->getName());
	}
}

void GEOObjects::mergeGeometries (std::vector<std::string> const & geo_names, std::string &merged_geo_name)
{
	std::vector<size_t> pnt_offsets(geo_names.size(), 0);

	// *** merge points
	std::vector<GEOLIB::Point*>* merged_points (new std::vector<GEOLIB::Point*>);
	for (size_t j(0); j<geo_names.size(); j++) {
		const std::vector<GEOLIB::Point*>* pnts (this->getPointVec(geo_names[j]));
		if (pnts) {
			size_t nPoints(0);
			// do not consider stations
			if (dynamic_cast<GEOLIB::Station*>((*pnts)[0]) == NULL) {
				nPoints = pnts->size();
				for (size_t k(0); k<nPoints; k++) {
					merged_points->push_back (new GEOLIB::Point (((*pnts)[k])->getCoords()));
				}
			}
			if (geo_names.size()-1 > j)
				pnt_offsets[j+1] = nPoints + pnt_offsets[j];
		}
	}

	this->addPointVec (merged_points, merged_geo_name);
	std::vector<size_t> const& id_map (this->getPointVecObj(merged_geo_name)->getIDMap ());

	// *** merge polylines
	std::vector<GEOLIB::Polyline*> *merged_polylines (new std::vector<GEOLIB::Polyline*>);
	for (size_t j(0); j<geo_names.size(); j++) {
		const std::vector<GEOLIB::Polyline*>* plys (this->getPolylineVec(geo_names[j]));
		if (plys) {
			for (size_t k(0); k<plys->size(); k++) {
				GEOLIB::Polyline* kth_ply_new(new GEOLIB::Polyline (*merged_points));
				GEOLIB::Polyline const*const kth_ply_old ((*plys)[k]);
				const size_t size_of_kth_ply (kth_ply_old->getNumberOfPoints());
				// copy point ids from old ply to new ply (considering the offset)
				for (size_t i(0); i<size_of_kth_ply; i++) {
					kth_ply_new->addPoint (id_map[pnt_offsets[j]+kth_ply_old->getPointID(i)]);
				}
				merged_polylines->push_back (kth_ply_new);
			}
		}
	}
	this->addPolylineVec (merged_polylines, merged_geo_name);


	// *** merge surfaces
	std::vector<GEOLIB::Surface*> *merged_sfcs (new std::vector<GEOLIB::Surface*>);
	for (size_t j(0); j<geo_names.size(); j++) {
		const std::vector<GEOLIB::Surface*>* sfcs (this->getSurfaceVec(geo_names[j]));
		if (sfcs) {
			for (size_t k(0); k<sfcs->size(); k++) {
				GEOLIB::Surface* kth_sfc_new(new GEOLIB::Surface (*merged_points));
				GEOLIB::Surface const*const kth_sfc_old ((*sfcs)[k]);
				const size_t size_of_kth_sfc (kth_sfc_old->getNTriangles());
				// copy point ids from old ply to new ply (considering the offset)
				for (size_t i(0); i<size_of_kth_sfc; i++) {
					const GEOLIB::Triangle* tri ((*kth_sfc_old)[i]);
					const size_t id0 (id_map[pnt_offsets[j]+(*tri)[0]]);
					const size_t id1 (id_map[pnt_offsets[j]+(*tri)[1]]);
					const size_t id2 (id_map[pnt_offsets[j]+(*tri)[2]]);
					kth_sfc_new->addTriangle (id0, id1, id2);
				}
				merged_sfcs->push_back (kth_sfc_new);
			}
		}
	}
	this->addSurfaceVec (merged_sfcs, merged_geo_name);
}


} // namespace
