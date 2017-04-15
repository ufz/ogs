/**
 * \file
 * \author Thomas Fischer / Karsten Rink
 * \date   2010-01-21
 * \brief  Implementation of the GEOObjects class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <fstream>

#include <logog/include/logog.hpp>

// GeoLib
#include "GEOObjects.h"

#include "Triangle.h"

namespace GeoLib
{
GEOObjects::GEOObjects() = default;

GEOObjects::~GEOObjects()
{
    // delete surfaces
    for (auto & _sfc_vec : _sfc_vecs)
        delete _sfc_vec;
    // delete polylines
    for (auto & _ply_vec : _ply_vecs)
        delete _ply_vec;
    // delete points
    for (auto & _pnt_vec : _pnt_vecs)
        delete _pnt_vec;
}

void GEOObjects::addPointVec(
    std::unique_ptr<std::vector<Point*>> points,
    std::string& name,
    std::map<std::string, std::size_t>* pnt_id_name_map,
    double eps)
{
    isUniquePointVecName(name);
    if (!points || points->empty()) {
        DBUG("GEOObjects::addPointVec(): Failed to create PointVec, because "
            "there aren't any points in the given vector.");
        return;
    }
    _pnt_vecs.push_back(new PointVec(name, std::move(points), pnt_id_name_map, PointVec::PointType::POINT, eps));
    _callbacks->addPointVec(name);
}

const std::vector<Point*>* GEOObjects::getPointVec(const std::string &name) const
{
    std::size_t const idx = this->exists(name);
    if (idx != std::numeric_limits<std::size_t>::max())
        return _pnt_vecs[idx]->getVector();

    DBUG("GEOObjects::getPointVec() - No entry found with name \"%s\".", name.c_str());
    return nullptr;
}

const PointVec* GEOObjects::getPointVecObj(const std::string &name) const
{
    std::size_t const idx = this->exists(name);
    if (idx != std::numeric_limits<std::size_t>::max())
        return _pnt_vecs[idx];

    DBUG("GEOObjects::getPointVecObj() - No entry found with name \"%s\".", name.c_str());
    return nullptr;
}

bool GEOObjects::removePointVec(std::string const& name)
{
    if (isPntVecUsed (name))
    {
        DBUG("GEOObjects::removePointVec() - There are still Polylines or Surfaces depending on these points.");
        return false;
    }

    for (std::vector<PointVec*>::iterator it(_pnt_vecs.begin());
         it != _pnt_vecs.end(); ++it)
        if ((*it)->getName().compare(name) == 0)
        {
            _callbacks->removePointVec(name);
            delete *it;
            _pnt_vecs.erase(it);
            return true;
        }
    DBUG("GEOObjects::removePointVec() - No entry found with name \"%s\".", name.c_str());
    return false;
}

void GEOObjects::addStationVec(std::unique_ptr<std::vector<Point*>> stations,
                               std::string& name)
{
    isUniquePointVecName(name);
    _pnt_vecs.push_back(new PointVec(name, std::move(stations), nullptr, PointVec::PointType::STATION));
    _callbacks->addStationVec(name);
}

const std::vector<GeoLib::Point*>* GEOObjects::getStationVec(
    const std::string& name) const
{
    for (std::vector<PointVec*>::const_iterator it(_pnt_vecs.begin());
         it != _pnt_vecs.end(); ++it) {
        if ((*it)->getName().compare(name) == 0 && (*it)->getType() == PointVec::PointType::STATION) {
            return (*it)->getVector();
        }
    }
    DBUG("GEOObjects::getStationVec() - No entry found with name \"%s\".", name.c_str());
    return nullptr;
}

void GEOObjects::addPolylineVec(std::unique_ptr<std::vector<Polyline*>> lines,
                                const std::string& name,
                                std::map<std::string, std::size_t>* ply_names)
{
    assert(lines);
    for (std::vector<Polyline*>::iterator it (lines->begin());
         it != lines->end(); )
    {
        if ((*it)->getNumberOfPoints() < 2)
        {
            std::vector<Polyline*>::iterator it_erase (it);
            it = lines->erase (it_erase);
        }
        else
            ++it;
    }

    if (lines->empty())
        return;

    _ply_vecs.push_back(new PolylineVec(name, std::move(lines), ply_names));
    _callbacks->addPolylineVec(name);
}

bool GEOObjects::appendPolylineVec(const std::vector<Polyline*>& polylines,
                                   const std::string& name)
{
    // search vector
    std::size_t idx (0);
    bool nfound (true);
    for (idx = 0; idx < _ply_vecs.size() && nfound; idx++)
        if ( (_ply_vecs[idx]->getName()).compare (name) == 0 )
            nfound = false;

    if (!nfound)
    {
        idx--;
        std::size_t n_plys (polylines.size());
        // append lines
        for (std::size_t k(0); k < n_plys; k++)
            _ply_vecs[idx]->push_back (polylines[k]);
        _callbacks->appendPolylineVec(name);
        return true;
    }
    else
        return false;
}

const std::vector<Polyline*>* GEOObjects::getPolylineVec(const std::string &name) const
{
    std::size_t size (_ply_vecs.size());
    for (std::size_t i = 0; i < size; i++)
        if (_ply_vecs[i]->getName().compare(name) == 0)
            return _ply_vecs[i]->getVector();

    DBUG("GEOObjects::getPolylineVec() - No entry found with name \"%s\".", name.c_str());
    return nullptr;
}

const PolylineVec* GEOObjects::getPolylineVecObj(const std::string &name) const
{
    std::size_t size (_ply_vecs.size());
    for (std::size_t i = 0; i < size; i++)
        if (_ply_vecs[i]->getName().compare(name) == 0)
            return _ply_vecs[i];

    DBUG("GEOObjects::getPolylineVecObj() - No entry found with name \"%s\".", name.c_str());
    return nullptr;
}

bool GEOObjects::removePolylineVec(std::string const& name)
{
    _callbacks->removePolylineVec(name);
    for (std::vector<PolylineVec*>::iterator it = _ply_vecs.begin();
         it != _ply_vecs.end(); ++it)
        if ((*it)->getName().compare(name) == 0)
        {
            delete *it;
            _ply_vecs.erase(it);
            return true;
        }

    DBUG("GEOObjects::removePolylineVec() - No entry found with name \"%s\".", name.c_str());
    return false;
}

void GEOObjects::addSurfaceVec(std::unique_ptr<std::vector<Surface*>> sfc,
                               const std::string& name,
                               std::map<std::string, std::size_t>* sfc_names)
{
    _sfc_vecs.push_back(new SurfaceVec(name, std::move(sfc), sfc_names));
    if (!sfc || !sfc->empty()) _callbacks->addSurfaceVec(name);
}

bool GEOObjects::appendSurfaceVec(const std::vector<Surface*>& surfaces,
                                  const std::string& name)
{
    // search vector
    std::size_t idx (0);
    bool nfound (true);
    for (idx = 0; idx < _sfc_vecs.size() && nfound; idx++)
        if ( (_sfc_vecs[idx]->getName()).compare (name) == 0 )
            nfound = false;

    if (!nfound)
    {
        idx--;
        std::size_t n_sfcs (surfaces.size());
        // append surfaces
        for (std::size_t k(0); k < n_sfcs; k++)
            _sfc_vecs[idx]->push_back (surfaces[k]);
        _callbacks->appendSurfaceVec(name);
        return true;
    }
    else
    {
        // the copy is needed because addSurfaceVec is passing it to SurfaceVec
        // ctor, which needs write access to the surface vector.
        auto sfc = std::unique_ptr<std::vector<GeoLib::Surface*>>{
            new std::vector<GeoLib::Surface*>};
        for (auto surface : surfaces)
            sfc->push_back(surface);
        addSurfaceVec(std::move(sfc), name);
        return false;
    }
}

const std::vector<Surface*>* GEOObjects::getSurfaceVec(const std::string &name) const
{
    std::size_t size (_sfc_vecs.size());
    for (std::size_t i = 0; i < size; i++)
        if (_sfc_vecs[i]->getName().compare(name) == 0)
            return _sfc_vecs[i]->getVector();
    DBUG("GEOObjects::getSurfaceVec() - No entry found with name \"%s\".", name.c_str());
    return nullptr;
}

bool GEOObjects::removeSurfaceVec(const std::string &name)
{
    _callbacks->removeSurfaceVec(name);
    for (std::vector<SurfaceVec*>::iterator it (_sfc_vecs.begin());
         it != _sfc_vecs.end(); ++it)
        if ((*it)->getName().compare (name) == 0)
        {
            delete *it;
            _sfc_vecs.erase (it);
            return true;
        }

    DBUG("GEOObjects::removeSurfaceVec() - No entry found with name \"%s\".", name.c_str());
    return false;
}

const SurfaceVec* GEOObjects::getSurfaceVecObj(const std::string &name) const
{
    std::size_t size (_sfc_vecs.size());
    for (std::size_t i = 0; i < size; i++)
        if (_sfc_vecs[i]->getName().compare(name) == 0)
            return _sfc_vecs[i];
    DBUG("GEOObjects::getSurfaceVecObj() - No entry found with name \"%s\".", name.c_str());
    return nullptr;
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
            cpName = cpName + "-" + std::to_string(count);

        for (auto & _pnt_vec : _pnt_vecs)
            if ( cpName.compare(_pnt_vec->getName()) == 0 )
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
         ++it)
    {
        if (((*it)->getName()).compare(name) == 0)
            return true;
    }
    for (std::vector<SurfaceVec*>::const_iterator it ( _sfc_vecs.begin()); it != _sfc_vecs.end();
         ++it)
    {
        if (((*it)->getName()).compare(name) == 0)
            return true;
    }

    return false;
}

void GEOObjects::getStationVectorNames(std::vector<std::string>& names) const
{
    for (std::vector<PointVec*>::const_iterator it(_pnt_vecs.begin()); it != _pnt_vecs.end();
         ++it)
        if ((*it)->getType() == PointVec::PointType::STATION)
            names.push_back((*it)->getName());
}

void GEOObjects::getGeometryNames (std::vector<std::string>& names) const
{
    names.clear ();
    for (std::vector<PointVec*>::const_iterator it(_pnt_vecs.begin()); it != _pnt_vecs.end();
         ++it)
        if ((*it)->getType() == PointVec::PointType::POINT)
            names.push_back((*it)->getName());
}

const std::string GEOObjects::getElementNameByID(const std::string &geometry_name, GeoLib::GEOTYPE type, std::size_t id) const
{
    std::string name("");
    switch (type)
    {
        case GeoLib::GEOTYPE::POINT:
            this->getPointVecObj(geometry_name)->getNameOfElementByID(id, name);
            break;
        case GeoLib::GEOTYPE::POLYLINE:
            this->getPolylineVecObj(geometry_name)->getNameOfElementByID(id, name);
            break;
        case GeoLib::GEOTYPE::SURFACE:
            this->getSurfaceVecObj(geometry_name)->getNameOfElementByID(id, name);
    }
    return name;
}

int GEOObjects::mergeGeometries (std::vector<std::string> const & geo_names,
                                 std::string &merged_geo_name)
{
    const std::size_t n_geo_names(geo_names.size());

    if (n_geo_names < 2)
        return 2;

    std::vector<std::size_t> pnt_offsets(n_geo_names, 0);

    if (! mergePoints(geo_names, merged_geo_name, pnt_offsets))
        return 1;

    mergePolylines(geo_names, merged_geo_name, pnt_offsets);

    mergeSurfaces(geo_names, merged_geo_name, pnt_offsets);

    return 0;
}

bool GEOObjects::mergePoints(std::vector<std::string> const & geo_names,
        std::string & merged_geo_name, std::vector<std::size_t> &pnt_offsets)
{
    const std::size_t n_geo_names(geo_names.size());

    auto merged_points = std::unique_ptr<std::vector<GeoLib::Point*>>(
        new std::vector<GeoLib::Point*>);
    std::map<std::string, std::size_t>* merged_pnt_names(new std::map<std::string, std::size_t>);

    for (std::size_t j(0); j < n_geo_names; ++j) {
        GeoLib::PointVec const*const pnt_vec(this->getPointVecObj(geo_names[j]));
        if (pnt_vec == nullptr)
            continue;
        const std::vector<GeoLib::Point*>* pnts(pnt_vec->getVector());
        if (pnts == nullptr) {
            return false;
        }

        // do not consider stations
        if (dynamic_cast<GeoLib::Station*>((*pnts)[0])) {
            continue;
        }

        std::size_t const n_pnts(pnts->size());
        for (std::size_t k(0); k < n_pnts; ++k) {
            merged_points->push_back(
                new GeoLib::Point(*(*pnts)[k], pnt_offsets[j]+k));
            std::string const& item_name(pnt_vec->getItemNameByID(k));
            if (! item_name.empty()) {
                merged_pnt_names->insert(
                    std::make_pair(item_name, pnt_offsets[j] + k));
            }
        }
        if (n_geo_names - 1 > j) {
            pnt_offsets[j + 1] = n_pnts + pnt_offsets[j];
        }
    }

    addPointVec (std::move(merged_points), merged_geo_name, merged_pnt_names, 1e-6);
    return true;
}

void GEOObjects::mergePolylines(std::vector<std::string> const & geo_names,
        std::string & merged_geo_name, std::vector<std::size_t> const& pnt_offsets)
{
    const std::size_t n_geo_names(geo_names.size());
    std::vector<std::size_t> ply_offsets(n_geo_names, 0);

    auto merged_polylines = std::unique_ptr<std::vector<GeoLib::Polyline*>>(
        new std::vector<GeoLib::Polyline*>);
    std::map<std::string, std::size_t>* merged_ply_names(new std::map<std::string, std::size_t>);

    std::vector<GeoLib::Point*> const* merged_points(this->getPointVecObj(merged_geo_name)->getVector());
    std::vector<std::size_t> const& id_map (this->getPointVecObj(merged_geo_name)->getIDMap ());

    for (std::size_t j(0); j < n_geo_names; j++) {
        const std::vector<GeoLib::Polyline*>* plys (this->getPolylineVec(geo_names[j]));
        if (plys) {
            std::string tmp_name;
            for (std::size_t k(0); k < plys->size(); k++) {
                GeoLib::Polyline* kth_ply_new(new GeoLib::Polyline (*merged_points));
                GeoLib::Polyline const* const kth_ply_old ((*plys)[k]);
                const std::size_t size_of_kth_ply (kth_ply_old->getNumberOfPoints());
                // copy point ids from old ply to new ply (considering the offset)
                for (std::size_t i(0); i < size_of_kth_ply; i++) {
                    kth_ply_new->addPoint (id_map[pnt_offsets[j] +
                                                  kth_ply_old->getPointID(i)]);
                }
                merged_polylines->push_back (kth_ply_new);
                if (this->getPolylineVecObj(geo_names[j])->getNameOfElementByID(k, tmp_name)) {
                    merged_ply_names->insert(std::pair<std::string, std::size_t>(tmp_name, ply_offsets[j] + k));
                }
            }
            if (n_geo_names - 1 > j) {
                ply_offsets[j + 1] = plys->size() + ply_offsets[j];
            }
        }
    }

    if (! merged_polylines->empty()) {
        this->addPolylineVec (std::move(merged_polylines), merged_geo_name, merged_ply_names);
    } else {
        delete merged_ply_names;
    }
}

void GEOObjects::mergeSurfaces(std::vector<std::string> const & geo_names,
        std::string & merged_geo_name, std::vector<std::size_t> const& pnt_offsets)
{
    std::vector<GeoLib::Point*> const* merged_points(this->getPointVecObj(merged_geo_name)->getVector());
    std::vector<std::size_t> const& id_map (this->getPointVecObj(merged_geo_name)->getIDMap ());

    const std::size_t n_geo_names(geo_names.size());
    std::vector<std::size_t> sfc_offsets(n_geo_names, 0);
    auto merged_sfcs = std::unique_ptr<std::vector<GeoLib::Surface*>>(
        new std::vector<GeoLib::Surface*>);
    std::map<std::string, std::size_t>* merged_sfc_names(new std::map<std::string, std::size_t>);
    for (std::size_t j(0); j < n_geo_names; j++) {
        const std::vector<GeoLib::Surface*>* sfcs (this->getSurfaceVec(geo_names[j]));
        if (sfcs) {
            std::string tmp_name;
            for (std::size_t k(0); k < sfcs->size(); k++) {
                GeoLib::Surface* kth_sfc_new(new GeoLib::Surface (*merged_points));
                GeoLib::Surface const* const kth_sfc_old ((*sfcs)[k]);
                const std::size_t size_of_kth_sfc (kth_sfc_old->getNumberOfTriangles());
                // clone surface elements using new ids
                for (std::size_t i(0); i < size_of_kth_sfc; i++) {
                    const GeoLib::Triangle* tri ((*kth_sfc_old)[i]);
                    const std::size_t id0 (id_map[pnt_offsets[j] + (*tri)[0]]);
                    const std::size_t id1 (id_map[pnt_offsets[j] + (*tri)[1]]);
                    const std::size_t id2 (id_map[pnt_offsets[j] + (*tri)[2]]);
                    kth_sfc_new->addTriangle (id0, id1, id2);
                }
                merged_sfcs->push_back (kth_sfc_new);

                if (this->getSurfaceVecObj(geo_names[j])->getNameOfElementByID(k, tmp_name)) {
                    merged_sfc_names->insert(std::pair<std::string, std::size_t>(tmp_name, sfc_offsets[j] + k));
                }
            }
            if (n_geo_names - 1 > j) {
                sfc_offsets[j + 1] = sfcs->size() + sfc_offsets[j];
            }
        }
    }
    if (! merged_sfcs->empty()) {
        this->addSurfaceVec (std::move(merged_sfcs), merged_geo_name, merged_sfc_names);
    } else {
        delete merged_sfc_names;
    }
}

const GeoLib::GeoObject* GEOObjects::getGeoObject(const std::string &geo_name,
                                                            GeoLib::GEOTYPE type,
                                                            const std::string &geo_obj_name) const
{
    GeoLib::GeoObject *geo_obj(nullptr);
    switch (type) {
    case GeoLib::GEOTYPE::POINT: {
        GeoLib::PointVec const* pnt_vec(getPointVecObj(geo_name));
        if (pnt_vec)
            geo_obj = const_cast<GeoLib::GeoObject*>(
                    dynamic_cast<GeoLib::GeoObject const*>(
                            pnt_vec->getElementByName(geo_obj_name)));
        break;
    }
    case GeoLib::GEOTYPE::POLYLINE: {
        GeoLib::PolylineVec const* ply_vec(getPolylineVecObj(geo_name));
        if (ply_vec)
            geo_obj = const_cast<GeoLib::GeoObject*>(
                    dynamic_cast<GeoLib::GeoObject const*>(
                            ply_vec->getElementByName(geo_obj_name)));
        break;
    }
    case GeoLib::GEOTYPE::SURFACE: {
        GeoLib::SurfaceVec const* sfc_vec(getSurfaceVecObj(geo_name));
        if (sfc_vec)
            geo_obj = const_cast<GeoLib::GeoObject*>(
                    dynamic_cast<GeoLib::GeoObject const*>(
                            sfc_vec->getElementByName(geo_obj_name)));
        break;
    }
    default:
        ERR("GEOObjects::getGeoObject(): geometric type not handled.")
        return nullptr;
    };

    if (!geo_obj) {
        DBUG("GEOObjects::getGeoObject(): Could not find %s \"%s\" in geometry.",
                GeoLib::convertGeoTypeToString(type).c_str(), geo_obj_name.c_str());
    }
    return geo_obj;
}

GeoLib::GeoObject const* GEOObjects::getGeoObject(
    const std::string &geo_name,
    const std::string &geo_obj_name) const
{
    GeoLib::GeoObject const* geo_obj(
        getGeoObject(geo_name, GeoLib::GEOTYPE::POINT, geo_obj_name)
    );

    if(!geo_obj)
        geo_obj = getGeoObject(geo_name, GeoLib::GEOTYPE::POLYLINE, geo_obj_name);

    if(!geo_obj)
        geo_obj = getGeoObject(geo_name, GeoLib::GEOTYPE::SURFACE, geo_obj_name);

    if (!geo_obj) {
        DBUG("GEOObjects::getGeoObject(): Could not find \"%s\" in geometry %s.",
            geo_obj_name.c_str(), geo_name.c_str());
    }
    return geo_obj;
}

std::size_t GEOObjects::exists(const std::string &geometry_name) const
{
    std::size_t const size (_pnt_vecs.size());
    for (std::size_t i = 0; i < size; i++)
        if (_pnt_vecs[i]->getName().compare(geometry_name) == 0)
            return i;

    // HACK for enabling conversion of files without loading the associated geometry
    if (size>0 && _pnt_vecs[0]->getName().compare("conversionTestRun#1")==0)
        return 1;

    return std::numeric_limits<std::size_t>::max();
}


} // namespace
