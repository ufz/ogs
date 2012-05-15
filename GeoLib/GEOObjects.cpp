/*
 * GEOObjects.cpp
 *
 *  Created on: Jan 21, 2010
 *      Author: TF / KR
 */

#include "GEOObjects.h"
#include "StringTools.h"

#include <fstream>

namespace GeoLib
{
GEOObjects::GEOObjects()
{
}

GEOObjects::~GEOObjects()
{
	// delete surfaces
	for (size_t k(0); k < _sfc_vecs.size(); k++)
		delete _sfc_vecs[k];
	// delete polylines
	for (size_t k(0); k < _ply_vecs.size(); k++)
		delete _ply_vecs[k];
	// delete points
	for (size_t k(0); k < _pnt_vecs.size(); k++)
		delete _pnt_vecs[k];
}

void GEOObjects::addPointVec(std::vector<Point*>* points,
                             std::string &name,
                             std::map<std::string, size_t>* pnt_id_name_map,
                             double eps)
{
	isUniquePointVecName(name);
	_pnt_vecs.push_back(new PointVec(name, points, pnt_id_name_map, PointVec::POINT, eps));
}

bool GEOObjects::appendPointVec(std::vector<Point*> const& new_points,
                                std::string const &name, std::vector<size_t>* ids)
{
	// search vector
	int idx = this->exists(name);

	if (idx>=0) {
		size_t n_new_pnts (new_points.size());
		// append points
		if (ids)
			for (size_t k(0); k < n_new_pnts; k++)
				ids->push_back (_pnt_vecs[idx]->push_back (new_points[k]));
		else
			for (size_t k(0); k < n_new_pnts; k++)
				_pnt_vecs[idx]->push_back (new_points[k]);

		return true;
	} else
		return false;
}

bool GEOObjects::appendPoint(Point* point, std::string const &name, size_t& id)
{
	// search vector
	int idx = this->exists(name);

	if (idx>=0) {
		const size_t size_previous (_pnt_vecs[idx]->size());
		id = _pnt_vecs[idx]->push_back (point);
		if (size_previous < _pnt_vecs[idx]->size()) {
			return true;
		} else {
			return false;
		}
	} else {
		return false;
	}
}

const std::vector<Point*>* GEOObjects::getPointVec(const std::string &name) const
{
	int idx = this->exists(name);
	if (idx>=0) return _pnt_vecs[idx]->getVector();
/*
	size_t size (_pnt_vecs.size());
	for (size_t i = 0; i < size; i++)
		if (_pnt_vecs[i]->getName().compare(name) == 0)
			return _pnt_vecs[i]->getVector();
*/
	std::cout << "GEOObjects::getPointVec() - No entry found with name \"" << name << "\"." << std::endl;
	return NULL;
}

const PointVec* GEOObjects::getPointVecObj(const std::string &name) const
{
	int idx = this->exists(name);
	if (idx>=0) return _pnt_vecs[idx];
/*
	size_t size (_pnt_vecs.size());
	for (size_t i = 0; i < size; i++)
		if (_pnt_vecs[i]->getName().compare(name) == 0)
			return _pnt_vecs[i];
*/
	std::cout << "GEOObjects::getPointVecObj() - No entry found with name \"" << name << "\"." << std::endl;
	return NULL;
}

bool GEOObjects::removePointVec(const std::string &name)
{
	if (isPntVecUsed (name))
	{
		std::cout <<
		"GEOObjects::removePointVec() - There are still Polylines or Surfaces depending on these points."
		          << std::endl;
		return false;
	}

	for (std::vector<PointVec*>::iterator it(_pnt_vecs.begin());
	     it != _pnt_vecs.end(); it++)
		if ((*it)->getName().compare(name) == 0)
		{
			delete *it;
			_pnt_vecs.erase(it);
			return true;
		}
	std::cout << "GEOObjects::removePointVec() - No entry found with name \"" << name << "." <<
	std::endl;
	return false;
}

void GEOObjects::addStationVec(std::vector<Point*>* stations, std::string &name)
{
	isUniquePointVecName(name);
	_pnt_vecs.push_back(new PointVec(name, stations, NULL, PointVec::STATION));
}

std::vector<Point*>* GEOObjects::filterStationVec(const std::string &name,
                                                  const std::vector<PropertyBounds> &bounds)
{
	for (std::vector<PointVec*>::iterator it(_pnt_vecs.begin());
	     it != _pnt_vecs.end(); it++)
		if ((*it)->getName().compare(name) == 0 && (*it)->getType()
		    == PointVec::STATION)
			return (*it)->filterStations(bounds);

	std::cout << "GEOObjects:: filterStations() - No entry found with name \""
	          << name << "." << std::endl;
	return NULL;
}

const std::vector<Point*>* GEOObjects::getStationVec(const std::string &name) const
{
	for (std::vector<PointVec*>::const_iterator it(_pnt_vecs.begin());
	     it != _pnt_vecs.end(); it++) {
		if ((*it)->getName().compare(name) == 0 && (*it)->getType() == PointVec::STATION) {
			return (*it)->getVector();
		}
	}
	std::cout << "GEOObjects::getStationVec() - No entry found with name \""
	          << name << "." << std::endl;
	return NULL;
}

void GEOObjects::addPolylineVec(std::vector<Polyline*>* lines,
                                const std::string &name, std::map<std::string, size_t>* ply_names)
{
	for (std::vector<Polyline*>::iterator it (lines->begin());
	     it != lines->end(); )
	{
		if ((*it)->getNumberOfPoints() < 2)
		{
			std::vector<Polyline*>::iterator it_erase (it);
			it = lines->erase (it_erase);
		}
		else
			it++;
	}

	if (lines->empty())
		return;

	_ply_vecs.push_back(new PolylineVec(name, lines, ply_names));
}

bool GEOObjects::appendPolylineVec(const std::vector<Polyline*> &polylines, const std::string &name)
{
	// search vector
	size_t idx (0);
	bool nfound (true);
	for (idx = 0; idx < _ply_vecs.size() && nfound; idx++)
		if ( (_ply_vecs[idx]->getName()).compare (name) == 0 )
			nfound = false;

	if (!nfound)
	{
		idx--;
		size_t n_plys (polylines.size());
		// append lines
		for (size_t k(0); k < n_plys; k++)
			_ply_vecs[idx]->push_back (polylines[k]);
		return true;
	}
	else
		return false;
}

const std::vector<Polyline*>* GEOObjects::getPolylineVec(const std::string &name) const
{
	size_t size (_ply_vecs.size());
	for (size_t i = 0; i < size; i++)
		if (_ply_vecs[i]->getName().compare(name) == 0)
			return _ply_vecs[i]->getVector();

#ifndef NDEBUG
	std::cout << "DEB: GEOObjects::getPolylineVec() - No entry found with name \"" << name <<
	"." << std::endl;
#endif
	return NULL;
}

const PolylineVec* GEOObjects::getPolylineVecObj(const std::string &name) const
{
	size_t size (_ply_vecs.size());
	for (size_t i = 0; i < size; i++)
		if (_ply_vecs[i]->getName().compare(name) == 0)
			return _ply_vecs[i];

#ifndef NDEBUG
	std::cout << "DEB: GEOObjects::getPolylineVecObj() - No entry found with name \"" << name <<
	"\"." << std::endl;
#endif
	return NULL;
}

bool GEOObjects::removePolylineVec(const std::string &name)
{
	for (std::vector<PolylineVec*>::iterator it = _ply_vecs.begin();
	     it != _ply_vecs.end(); ++it)
		if ((*it)->getName().compare(name) == 0)
		{
			delete *it;
			_ply_vecs.erase(it);
			return true;
		}

#ifndef NDEBUG
	std::cout << "GEOObjects::removePolylineVec() - No entry found with name \""
	          << name << "\"." << std::endl;
#endif
	return false;
}

void GEOObjects::addSurfaceVec(std::vector<Surface*>* sfc, const std::string &name,
                               std::map<std::string, size_t>* sfc_names)
{
	_sfc_vecs.push_back(new SurfaceVec(name, sfc, sfc_names));
}

bool GEOObjects::appendSurfaceVec(const std::vector<Surface*> &surfaces, const std::string &name)
{
	// search vector
	size_t idx (0);
	bool nfound (true);
	for (idx = 0; idx < _sfc_vecs.size() && nfound; idx++)
		if ( (_sfc_vecs[idx]->getName()).compare (name) == 0 )
			nfound = false;

	if (!nfound)
	{
		idx--;
		size_t n_sfcs (surfaces.size());
		// append surfaces
		for (size_t k(0); k < n_sfcs; k++)
			_sfc_vecs[idx]->push_back (surfaces[k]);
		return true;
	}
	else
		return false;
}

const std::vector<Surface*>* GEOObjects::getSurfaceVec(const std::string &name) const
{
	size_t size (_sfc_vecs.size());
	for (size_t i = 0; i < size; i++)
		if (_sfc_vecs[i]->getName().compare(name) == 0)
			return _sfc_vecs[i]->getVector();
	std::cout << "GEOObjects::getSurfaceVec() - No entry found with name \"" << name << "." <<
	std::endl;
	return NULL;
}

bool GEOObjects::removeSurfaceVec(const std::string &name)
{
	for (std::vector<SurfaceVec*>::iterator it (_sfc_vecs.begin());
	     it != _sfc_vecs.end(); it++)
		if ((*it)->getName().compare (name) == 0)
		{
			delete *it;
			_sfc_vecs.erase (it);
			return true;
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
	for (size_t i = 0; i < size; i++)
		if (_sfc_vecs[i]->getName().compare(name) == 0)
			return _sfc_vecs[i];
	std::cout << "GEOObjects::getSurfaceVecObj() - No entry found with name \"" << name <<
	"\"." << std::endl;
	return NULL;
}

bool GEOObjects::isUniquePointVecName(std::string &name)
{
	int count = 0;
	bool isUnique = false;
	std::string cpName;

	while (!isUnique)
	{
		isUnique = true;
		cpName = name;

		count++;
		// If the original name already exists we start to add numbers to name for
		// as long as it takes to make the name unique.
		if (count > 1)
			cpName = cpName + "-" + number2str(count);

		for (size_t i = 0; i < _pnt_vecs.size(); i++)
			if ( cpName.compare(_pnt_vecs[i]->getName()) == 0 )
				isUnique = false;
	}

	// At this point cpName is a unique name and isUnique is true.
	// If cpName is not the original name, "name" is changed and isUnique is set to false,
	// indicating that a vector with the original name already exists.
	if (count > 1)
	{
		isUnique = false;
		name = cpName;
	}
	return isUnique;
}

bool GEOObjects::isPntVecUsed (const std::string &name) const
{
	// search dependent data structures (Polyline)
	for (std::vector<PolylineVec*>::const_iterator it ( _ply_vecs.begin()); it != _ply_vecs.end();
	     it++)
	{
		std::string a = (*it)->getName();
		if (((*it)->getName()).compare(name) == 0)
			return true;
	}
	for (std::vector<SurfaceVec*>::const_iterator it ( _sfc_vecs.begin()); it != _sfc_vecs.end();
	     it++)
	{
		std::string a = (*it)->getName();
		if (((*it)->getName()).compare(name) == 0)
			return true;
	}

	return false;
}

void GEOObjects::getStationVectorNames(std::vector<std::string>& names) const
{
	for (std::vector<PointVec*>::const_iterator it(_pnt_vecs.begin()); it != _pnt_vecs.end();
	     it++)
		if ((*it)->getType() == PointVec::STATION)
			names.push_back((*it)->getName());
}

void GEOObjects::getGeometryNames (std::vector<std::string>& names) const
{
	names.clear ();
	for (std::vector<PointVec*>::const_iterator it(_pnt_vecs.begin()); it != _pnt_vecs.end();
	     it++)
		if ((*it)->getType() == PointVec::POINT)
			names.push_back((*it)->getName());
}

const std::string GEOObjects::getElementNameByID(const std::string &geometry_name, GeoLib::GEOTYPE type, size_t id) const
{
	std::string name("");
	switch (type)
	{
		case GeoLib::POINT:
			this->getPointVecObj(geometry_name)->getNameOfElementByID(id, name);
			break;
		case GeoLib::POLYLINE:
			this->getPolylineVecObj(geometry_name)->getNameOfElementByID(id, name);
			break;
		case GeoLib::SURFACE:
			this->getSurfaceVecObj(geometry_name)->getNameOfElementByID(id, name);
			break;
		default:
			std::cout << "No valid GEOTYPE given." << std::endl;
	}
	return name;
}

void GEOObjects::mergeGeometries (std::vector<std::string> const & geo_names,
                                  std::string &merged_geo_name)
{
	const size_t n_geo_names(geo_names.size());
	std::vector<size_t> pnt_offsets(n_geo_names, 0);

	// *** merge points
	std::vector<GeoLib::Point*>* merged_points (new std::vector<GeoLib::Point*>);
	for (size_t j(0); j < n_geo_names; j++) {
		const std::vector<GeoLib::Point*>* pnts (this->getPointVec(geo_names[j]));
		if (pnts) {
			size_t n_pnts(0);
			// do not consider stations
			if (dynamic_cast<GeoLib::Station*>((*pnts)[0]) == NULL) {
				n_pnts = pnts->size();
				for (size_t k(0); k < n_pnts; k++)
					merged_points->push_back (new GeoLib::Point (((*pnts)[k])->getCoords()));
			}
			if (n_geo_names - 1 > j) {
				pnt_offsets[j + 1] = n_pnts + pnt_offsets[j];
			}
		}
	}
	addPointVec (merged_points, merged_geo_name, NULL, 1e-6);
	std::vector<size_t> const& id_map (this->getPointVecObj(merged_geo_name)->getIDMap ());

	// *** merge polylines
	std::vector<GeoLib::Polyline*>* merged_polylines (new std::vector<GeoLib::Polyline*>);
	for (size_t j(0); j < n_geo_names; j++) {
		const std::vector<GeoLib::Polyline*>* plys (this->getPolylineVec(geo_names[j]));
		if (plys) {
			for (size_t k(0); k < plys->size(); k++) {
				GeoLib::Polyline* kth_ply_new(new GeoLib::Polyline (*merged_points));
				GeoLib::Polyline const* const kth_ply_old ((*plys)[k]);
				const size_t size_of_kth_ply (kth_ply_old->getNumberOfPoints());
				// copy point ids from old ply to new ply (considering the offset)
				for (size_t i(0); i < size_of_kth_ply; i++) {
					kth_ply_new->addPoint (id_map[pnt_offsets[j] +
					                              kth_ply_old->getPointID(i)]);
				}
				merged_polylines->push_back (kth_ply_new);
			}
		}
	}
	this->addPolylineVec (merged_polylines, merged_geo_name);

	// *** merge surfaces
	std::vector<GeoLib::Surface*>* merged_sfcs (new std::vector<GeoLib::Surface*>);
	for (size_t j(0); j < n_geo_names; j++) {
		const std::vector<GeoLib::Surface*>* sfcs (this->getSurfaceVec(geo_names[j]));
		if (sfcs) {
			for (size_t k(0); k < sfcs->size(); k++) {
				GeoLib::Surface* kth_sfc_new(new GeoLib::Surface (*merged_points));
				GeoLib::Surface const* const kth_sfc_old ((*sfcs)[k]);
				const size_t size_of_kth_sfc (kth_sfc_old->getNTriangles());
				// copy point ids from old ply to new ply (considering the offset)
				for (size_t i(0); i < size_of_kth_sfc; i++) {
					const GeoLib::Triangle* tri ((*kth_sfc_old)[i]);
					const size_t id0 (id_map[pnt_offsets[j] + (*tri)[0]]);
					const size_t id1 (id_map[pnt_offsets[j] + (*tri)[1]]);
					const size_t id2 (id_map[pnt_offsets[j] + (*tri)[2]]);
					kth_sfc_new->addTriangle (id0, id1, id2);
				}
				merged_sfcs->push_back (kth_sfc_new);
			}
		}
	}
	this->addSurfaceVec (merged_sfcs, merged_geo_name);
}

const GeoLib::GeoObject* GEOObjects::getGEOObject(const std::string &geo_name,
                                                            GeoLib::GEOTYPE type,
                                                            const std::string &obj_name) const
{
	if (type == GeoLib::POINT)
		return this->getPointVecObj(geo_name)->getElementByName(obj_name);
	else if (type == GeoLib::POLYLINE)
		return this->getPolylineVecObj(geo_name)->getElementByName(obj_name);
	else if (type == GeoLib::SURFACE)
		return this->getSurfaceVecObj(geo_name)->getElementByName(obj_name);
	return NULL;
}

int GEOObjects::exists(const std::string &geometry_name) const
{
	size_t size (_pnt_vecs.size());
	for (size_t i = 0; i < size; i++)
		if (_pnt_vecs[i]->getName().compare(geometry_name) == 0)
			return i;

	// HACK for enabling conversion of files without loading the associated geometry
	if (size>0 && _pnt_vecs[0]->getName().compare("conversionTestRun#1")==0)	
		return 1;

	return -1;
}


} // namespace
