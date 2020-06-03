/**
 * \file
 * \author Thomas Fischer / Karsten Rink
 * \date   2010-01-21
 * \brief  Implementation of the GEOObjects class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GEOObjects.h"

#include <fstream>
#include "BaseLib/Logging.h"

#include "Triangle.h"

namespace GeoLib
{
GEOObjects::GEOObjects() = default;

GEOObjects::~GEOObjects()
{
    // delete surfaces
    for (auto& surface : sfc_vecs_)
    {
        delete surface;
    }
    // delete polylines
    for (auto& polyline : ply_vecs_)
    {
        delete polyline;
    }
    // delete points
    for (auto& point : pnt_vecs_)
    {
        delete point;
    }
}

void GEOObjects::addPointVec(std::unique_ptr<std::vector<Point*>> points,
                             std::string& name,
                             std::unique_ptr<std::map<std::string, std::size_t>>
                                 pnt_id_name_map,
                             double eps)
{
    isUniquePointVecName(name);
    if (!points || points->empty()) {
        DBUG("GEOObjects::addPointVec(): Failed to create PointVec, because "
            "there aren't any points in the given vector.");
        return;
    }
    pnt_vecs_.push_back(new PointVec(name, std::move(points),
                                     std::move(pnt_id_name_map),
                                     PointVec::PointType::POINT, eps));
    callbacks_->addPointVec(name);
}

const std::vector<Point*>* GEOObjects::getPointVec(const std::string &name) const
{
    std::size_t const idx = this->exists(name);
    if (idx != std::numeric_limits<std::size_t>::max())
    {
        return pnt_vecs_[idx]->getVector();
    }

    DBUG("GEOObjects::getPointVec() - No entry found with name '{:s}'.", name);
    return nullptr;
}

const PointVec* GEOObjects::getPointVecObj(const std::string &name) const
{
    std::size_t const idx = this->exists(name);
    if (idx != std::numeric_limits<std::size_t>::max())
    {
        return pnt_vecs_[idx];
    }

    DBUG("GEOObjects::getPointVecObj() - No entry found with name '{:s}'.",
         name);
    return nullptr;
}

bool GEOObjects::removePointVec(std::string const& name)
{
    if (isPntVecUsed (name))
    {
        DBUG("GEOObjects::removePointVec() - There are still Polylines or Surfaces depending on these points.");
        return false;
    }

    for (auto it(pnt_vecs_.begin()); it != pnt_vecs_.end(); ++it)
    {
        if ((*it)->getName() == name)
        {
            callbacks_->removePointVec(name);
            delete *it;
            pnt_vecs_.erase(it);
            return true;
        }
    }
    DBUG("GEOObjects::removePointVec() - No entry found with name '{:s}'.",
         name);
    return false;
}

void GEOObjects::addStationVec(std::unique_ptr<std::vector<Point*>> stations,
                               std::string& name)
{
    isUniquePointVecName(name);
    pnt_vecs_.push_back(new PointVec(name, std::move(stations), nullptr, PointVec::PointType::STATION));
    callbacks_->addStationVec(name);
}

const std::vector<GeoLib::Point*>* GEOObjects::getStationVec(
    const std::string& name) const
{
    auto const it = std::find_if(
        begin(pnt_vecs_), end(pnt_vecs_), [&name](PointVec const* const p) {
            return p->getName() == name &&
                   p->getType() == PointVec::PointType::STATION;
        });
    if (it != end(pnt_vecs_))
    {
        return (*it)->getVector();
    }
    DBUG("GEOObjects::getStationVec() - No entry found with name '{:s}'.",
         name);
    return nullptr;
}

void GEOObjects::addPolylineVec(
    std::unique_ptr<std::vector<Polyline*>> lines,
    const std::string& name,
    std::unique_ptr<std::map<std::string, std::size_t>>
        ply_names)
{
    assert(lines);
    for (auto it(lines->begin()); it != lines->end();)
    {
        if ((*it)->getNumberOfPoints() < 2)
        {
            auto it_erase(it);
            it = lines->erase (it_erase);
        }
        else
        {
            ++it;
        }
    }

    if (lines->empty())
    {
        return;
    }

    ply_vecs_.push_back(
        new PolylineVec(name, std::move(lines), std::move(ply_names)));
    callbacks_->addPolylineVec(name);
}

bool GEOObjects::appendPolylineVec(const std::vector<Polyline*>& polylines,
                                   const std::string& name)
{
    // search vector
    std::size_t idx (0);
    bool nfound (true);
    for (idx = 0; idx < ply_vecs_.size() && nfound; idx++)
    {
        if (ply_vecs_[idx]->getName() == name)
        {
            nfound = false;
        }
    }

    if (!nfound)
    {
        idx--;
        std::size_t n_plys (polylines.size());
        // append lines
        for (std::size_t k(0); k < n_plys; k++)
        {
            ply_vecs_[idx]->push_back(polylines[k]);
        }
        callbacks_->appendPolylineVec(name);
        return true;
    }

    return false;
}

const std::vector<Polyline*>* GEOObjects::getPolylineVec(const std::string &name) const
{
    std::size_t size (ply_vecs_.size());
    for (std::size_t i = 0; i < size; i++)
    {
        if (ply_vecs_[i]->getName() == name)
        {
            return ply_vecs_[i]->getVector();
        }
    }

    DBUG("GEOObjects::getPolylineVec() - No entry found with name '{:s}'.",
         name);
    return nullptr;
}

const PolylineVec* GEOObjects::getPolylineVecObj(const std::string &name) const
{
    std::size_t size (ply_vecs_.size());
    for (std::size_t i = 0; i < size; i++)
    {
        if (ply_vecs_[i]->getName() == name)
        {
            return ply_vecs_[i];
        }
    }

    DBUG("GEOObjects::getPolylineVecObj() - No entry found with name '{:s}'.",
         name);
    return nullptr;
}

bool GEOObjects::removePolylineVec(std::string const& name)
{
    callbacks_->removePolylineVec(name);
    for (auto it = ply_vecs_.begin(); it != ply_vecs_.end(); ++it)
    {
        if ((*it)->getName() == name)
        {
            delete *it;
            ply_vecs_.erase(it);
            return true;
        }
    }

    DBUG("GEOObjects::removePolylineVec() - No entry found with name '{:s}'.",
         name);
    return false;
}

void GEOObjects::addSurfaceVec(
    std::unique_ptr<std::vector<Surface*>> sfc,
    const std::string& name,
    std::unique_ptr<std::map<std::string, std::size_t>>
        sfc_names)
{
    sfc_vecs_.push_back(
        new SurfaceVec(name, std::move(sfc), std::move(sfc_names)));
    callbacks_->addSurfaceVec(name);
}

bool GEOObjects::appendSurfaceVec(const std::vector<Surface*>& surfaces,
                                  const std::string& name)
{
    // search vector
    std::size_t idx (0);
    bool nfound (true);
    for (idx = 0; idx < sfc_vecs_.size() && nfound; idx++)
    {
        if (sfc_vecs_[idx]->getName() == name)
        {
            nfound = false;
        }
    }

    if (!nfound)
    {
        idx--;
        std::size_t n_sfcs (surfaces.size());
        // append surfaces
        for (std::size_t k(0); k < n_sfcs; k++)
        {
            sfc_vecs_[idx]->push_back(surfaces[k]);
        }
        callbacks_->appendSurfaceVec(name);
        return true;
    }

    // the copy is needed because addSurfaceVec is passing it to SurfaceVec
    // ctor, which needs write access to the surface vector.
    auto sfc = std::make_unique<std::vector<GeoLib::Surface*>>(begin(surfaces),
                                                               end(surfaces));
    addSurfaceVec(std::move(sfc), name);
    return false;
}

const std::vector<Surface*>* GEOObjects::getSurfaceVec(const std::string &name) const
{
    std::size_t size (sfc_vecs_.size());
    for (std::size_t i = 0; i < size; i++)
    {
        if (sfc_vecs_[i]->getName() == name)
        {
            return sfc_vecs_[i]->getVector();
        }
    }
    DBUG("GEOObjects::getSurfaceVec() - No entry found with name '{:s}'.",
         name);
    return nullptr;
}

bool GEOObjects::removeSurfaceVec(const std::string &name)
{
    callbacks_->removeSurfaceVec(name);
    for (auto it(sfc_vecs_.begin()); it != sfc_vecs_.end(); ++it)
    {
        if ((*it)->getName() == name)
        {
            delete *it;
            sfc_vecs_.erase (it);
            return true;
        }
    }

    DBUG("GEOObjects::removeSurfaceVec() - No entry found with name '{:s}'.",
         name);
    return false;
}

const SurfaceVec* GEOObjects::getSurfaceVecObj(const std::string &name) const
{
    std::size_t size (sfc_vecs_.size());
    for (std::size_t i = 0; i < size; i++)
    {
        if (sfc_vecs_[i]->getName() == name)
        {
            return sfc_vecs_[i];
        }
    }
    DBUG("GEOObjects::getSurfaceVecObj() - No entry found with name '{:s}'.",
         name);
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
        {
            cpName = cpName + "-" + std::to_string(count);
        }

        for (auto& point : pnt_vecs_)
        {
            if (cpName == point->getName())
            {
                isUnique = false;
            }
        }
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
    for (auto polyline : ply_vecs_)
    {
        if (polyline->getName() == name)
        {
            return true;
        }
    }
    for (auto surface : sfc_vecs_)
    {
        if (surface->getName() == name)
        {
            return true;
        }
    }

    return false;
}

void GEOObjects::getStationVectorNames(std::vector<std::string>& names) const
{
    for (auto point : pnt_vecs_)
    {
        if (point->getType() == PointVec::PointType::STATION)
        {
            names.push_back(point->getName());
        }
    }
}

void GEOObjects::getGeometryNames (std::vector<std::string>& names) const
{
    names.clear ();
    for (auto point : pnt_vecs_)
    {
        if (point->getType() == PointVec::PointType::POINT)
        {
            names.push_back(point->getName());
        }
    }
}

std::string GEOObjects::getElementNameByID(const std::string& geometry_name,
                                           GeoLib::GEOTYPE type,
                                           std::size_t id) const
{
    std::string name;
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
    {
        return 2;
    }

    std::vector<std::size_t> pnt_offsets(n_geo_names, 0);

    if (!mergePoints(geo_names, merged_geo_name, pnt_offsets))
    {
        return 1;
    }

    mergePolylines(geo_names, merged_geo_name, pnt_offsets);

    mergeSurfaces(geo_names, merged_geo_name, pnt_offsets);

    return 0;
}

bool GEOObjects::mergePoints(std::vector<std::string> const & geo_names,
        std::string & merged_geo_name, std::vector<std::size_t> &pnt_offsets)
{
    const std::size_t n_geo_names(geo_names.size());

    auto merged_points = std::make_unique<std::vector<GeoLib::Point*>>();
    auto merged_pnt_names =
        std::make_unique<std::map<std::string, std::size_t>>();

    for (std::size_t j(0); j < n_geo_names; ++j) {
        GeoLib::PointVec const*const pnt_vec(this->getPointVecObj(geo_names[j]));
        if (pnt_vec == nullptr)
        {
            continue;
        }
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

    addPointVec(std::move(merged_points), merged_geo_name,
                std::move(merged_pnt_names), 1e-6);
    return true;
}

void GEOObjects::mergePolylines(std::vector<std::string> const& geo_names,
                                std::string const& merged_geo_name,
                                std::vector<std::size_t> const& pnt_offsets)
{
    const std::size_t n_geo_names(geo_names.size());
    std::vector<std::size_t> ply_offsets(n_geo_names, 0);

    auto merged_polylines = std::make_unique<std::vector<GeoLib::Polyline*>>();
    auto merged_ply_names =
        std::make_unique<std::map<std::string, std::size_t>>();

    std::vector<GeoLib::Point*> const* merged_points(this->getPointVecObj(merged_geo_name)->getVector());
    std::vector<std::size_t> const& id_map (this->getPointVecObj(merged_geo_name)->getIDMap ());

    for (std::size_t j(0); j < n_geo_names; j++) {
        const std::vector<GeoLib::Polyline*>* plys (this->getPolylineVec(geo_names[j]));
        if (plys) {
            std::string tmp_name;
            for (std::size_t k(0); k < plys->size(); k++) {
                auto* kth_ply_new(new GeoLib::Polyline(*merged_points));
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
        this->addPolylineVec(std::move(merged_polylines), merged_geo_name,
                             std::move(merged_ply_names));
    }
}

void GEOObjects::mergeSurfaces(std::vector<std::string> const& geo_names,
                               std::string const& merged_geo_name,
                               std::vector<std::size_t> const& pnt_offsets)
{
    std::vector<GeoLib::Point*> const* merged_points(this->getPointVecObj(merged_geo_name)->getVector());
    std::vector<std::size_t> const& id_map (this->getPointVecObj(merged_geo_name)->getIDMap ());

    const std::size_t n_geo_names(geo_names.size());
    std::vector<std::size_t> sfc_offsets(n_geo_names, 0);
    auto merged_sfcs = std::make_unique<std::vector<GeoLib::Surface*>>();
    std::unique_ptr<std::map<std::string, std::size_t>> merged_sfc_names;
    for (std::size_t j(0); j < n_geo_names; j++) {
        const std::vector<GeoLib::Surface*>* sfcs (this->getSurfaceVec(geo_names[j]));
        if (sfcs) {
            std::string tmp_name;
            for (std::size_t k(0); k < sfcs->size(); k++) {
                auto* kth_sfc_new(new GeoLib::Surface(*merged_points));
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
        this->addSurfaceVec(std::move(merged_sfcs), merged_geo_name,
                            std::move(merged_sfc_names));
    }
}

void GEOObjects::renameGeometry(std::string const& old_name,
                                std::string const& new_name)
{
    callbacks_->renameGeometry(old_name, new_name);
    for (auto* pnt_vec : pnt_vecs_)
    {
        if (pnt_vec->getName() == old_name)
        {
            pnt_vec->setName(new_name);
            break;
        }
    }
    for (auto* ply_vec : ply_vecs_)
    {
        if (ply_vec->getName() == old_name)
        {
            ply_vec->setName(new_name);
            break;
        }
    }
    for (auto* sfc_vec : sfc_vecs_)
    {
        if (sfc_vec->getName() == old_name)
        {
            sfc_vec->setName(new_name);
            break;
        }
    }
}

void GEOObjects::markUnusedPoints(std::string const& geo_name,
                                  std::vector<bool>& transfer_pnts) const
{
    GeoLib::PolylineVec const* const ply_obj(getPolylineVecObj(geo_name));
    if (ply_obj)
    {
        std::vector<GeoLib::Polyline*> const& lines(*ply_obj->getVector());
        for (auto* line : lines)
        {
            std::size_t const n_pnts(line->getNumberOfPoints());
            for (std::size_t i = 0; i < n_pnts; ++i)
            {
                transfer_pnts[line->getPointID(i)] = false;
            }
        }
    }

    GeoLib::SurfaceVec const* const sfc_obj(getSurfaceVecObj(geo_name));
    if (sfc_obj)
    {
        std::vector<GeoLib::Surface*> const& surfaces = *sfc_obj->getVector();
        for (auto* sfc : surfaces)
        {
            std::size_t const n_tri(sfc->getNumberOfTriangles());
            for (std::size_t i = 0; i < n_tri; ++i)
            {
                GeoLib::Triangle const& t = *(*sfc)[i];
                transfer_pnts[t[0]] = false;
                transfer_pnts[t[1]] = false;
                transfer_pnts[t[2]] = false;
            }
        }
    }
}

int GEOObjects::geoPointsToStations(std::string const& geo_name,
                                    std::string& stn_name,
                                    bool const only_unused_pnts)
{
    GeoLib::PointVec const* const pnt_obj(getPointVecObj(geo_name));
    if (pnt_obj == nullptr)
    {
        ERR("Point vector {:s} not found.", geo_name);
        return -1;
    }
    std::vector<GeoLib::Point*> const& pnts = *pnt_obj->getVector();
    if (pnts.empty())
    {
        ERR("Point vector {:s} is empty.", geo_name);
        return -1;
    }
    std::size_t const n_pnts(pnts.size());
    std::vector<bool> transfer_pnts(n_pnts, true);
    if (only_unused_pnts)
    {
        markUnusedPoints(geo_name, transfer_pnts);
    }

    auto stations = std::make_unique<std::vector<GeoLib::Point*>>();
    for (std::size_t i = 0; i < n_pnts; ++i)
    {
        if (!transfer_pnts[i])
        {
            continue;
        }
        std::string name = pnt_obj->getItemNameByID(i);
        if (name.empty())
        {
            name = "Station " + std::to_string(i);
        }
        stations->push_back(new GeoLib::Station(pnts[i], name));
    }
    if (!stations->empty())
    {
        addStationVec(std::move(stations), stn_name);
    }
    else
    {
        WARN("No points found to convert.");
        return 1;
    }
    return 0;
}

const GeoLib::GeoObject* GEOObjects::getGeoObject(
    const std::string& geo_name,
    GeoLib::GEOTYPE type,
    const std::string& geo_obj_name) const
{
    GeoLib::GeoObject *geo_obj(nullptr);
    switch (type) {
    case GeoLib::GEOTYPE::POINT: {
        GeoLib::PointVec const* pnt_vec(getPointVecObj(geo_name));
        if (pnt_vec)
        {
            geo_obj = const_cast<GeoLib::GeoObject*>(
                dynamic_cast<GeoLib::GeoObject const*>(
                    pnt_vec->getElementByName(geo_obj_name)));
        }
        break;
    }
    case GeoLib::GEOTYPE::POLYLINE: {
        GeoLib::PolylineVec const* ply_vec(getPolylineVecObj(geo_name));
        if (ply_vec)
        {
            geo_obj = const_cast<GeoLib::GeoObject*>(
                dynamic_cast<GeoLib::GeoObject const*>(
                    ply_vec->getElementByName(geo_obj_name)));
        }
        break;
    }
    case GeoLib::GEOTYPE::SURFACE: {
        GeoLib::SurfaceVec const* sfc_vec(getSurfaceVecObj(geo_name));
        if (sfc_vec)
        {
            geo_obj = const_cast<GeoLib::GeoObject*>(
                dynamic_cast<GeoLib::GeoObject const*>(
                    sfc_vec->getElementByName(geo_obj_name)));
        }
        break;
    }
    default:
        ERR("GEOObjects::getGeoObject(): geometric type not handled.");
        return nullptr;
    };

    if (!geo_obj) {
        DBUG(
            "GEOObjects::getGeoObject(): Could not find {:s} '{:s}' in "
            "geometry.",
            GeoLib::convertGeoTypeToString(type),
            geo_obj_name);
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

    if (!geo_obj)
    {
        geo_obj =
            getGeoObject(geo_name, GeoLib::GEOTYPE::POLYLINE, geo_obj_name);
    }

    if (!geo_obj)
    {
        geo_obj =
            getGeoObject(geo_name, GeoLib::GEOTYPE::SURFACE, geo_obj_name);
    }

    if (!geo_obj) {
        DBUG(
            "GEOObjects::getGeoObject(): Could not find '{:s}' in geometry "
            "{:s}.",
            geo_obj_name, geo_name);
    }
    return geo_obj;
}

std::size_t GEOObjects::exists(const std::string &geometry_name) const
{
    std::size_t const size (pnt_vecs_.size());
    for (std::size_t i = 0; i < size; i++)
    {
        if (pnt_vecs_[i]->getName() == geometry_name)
        {
            return i;
        }
    }

    // HACK for enabling conversion of files without loading the associated geometry
    if (size > 0 && pnt_vecs_[0]->getName() == "conversionTestRun#1")
    {
        return 1;
    }

    return std::numeric_limits<std::size_t>::max();
}

}  // namespace GeoLib
