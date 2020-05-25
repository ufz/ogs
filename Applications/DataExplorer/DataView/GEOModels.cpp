/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-09
 * \brief  Implementation of the GEOModels class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "GEOModels.h"

#include "BaseLib/Logging.h"

#include "Applications/FileIO/Legacy/createSurface.h"
#include "GeoTreeModel.h"
#include "StationTreeModel.h"

#include "OGSError.h"

#include "GeoLib/Triangle.h"

GEOModels::GEOModels(GeoLib::GEOObjects& geo_objects, QObject* parent /*= 0*/)
    : QObject(parent), geo_objects_(geo_objects)
{
    geo_objects_.callbacks_ = std::make_unique<GEOModelsCallbacks>(*this);

    geoModel_ = new GeoTreeModel();
    stationModel_ = new StationTreeModel();
}

GEOModels::~GEOModels()
{
    delete stationModel_;
    delete geoModel_;
}

void GEOModels::updateGeometry(const std::string &geo_name)
{
    GeoLib::PointVec* points (geo_objects_.getPointVecObj(geo_name));
    GeoLib::PolylineVec* lines (geo_objects_.getPolylineVecObj(geo_name));
    GeoLib::SurfaceVec* surfaces (geo_objects_.getSurfaceVecObj(geo_name));

    if (points)
    {
        emit geoDataRemoved(geoModel_, geo_name, GeoLib::GEOTYPE::POINT);
        this->geoModel_->removeGeoList(geo_name, GeoLib::GEOTYPE::POINT);
        geoModel_->addPointList(QString::fromStdString(geo_name), *points);
        emit geoDataAdded(geoModel_, geo_name, GeoLib::GEOTYPE::POINT);

        if (lines)
        {
            emit geoDataRemoved(geoModel_, geo_name, GeoLib::GEOTYPE::POLYLINE);
            this->geoModel_->removeGeoList(geo_name, GeoLib::GEOTYPE::POLYLINE);
            geoModel_->addPolylineList(QString::fromStdString(geo_name), *lines);
            emit geoDataAdded(geoModel_, geo_name, GeoLib::GEOTYPE::POLYLINE);
        }

        if (surfaces)
        {
            emit geoDataRemoved(geoModel_, geo_name, GeoLib::GEOTYPE::SURFACE);
            this->geoModel_->removeGeoList(geo_name, GeoLib::GEOTYPE::SURFACE);
            geoModel_->addSurfaceList(QString::fromStdString(geo_name), *surfaces);
            emit geoDataAdded(geoModel_, geo_name, GeoLib::GEOTYPE::SURFACE);
        }
    }
    else
        ERR("GEOModels::updateGeometry() - Geometry '{:s}' not found.",
            geo_name);
}

void GEOModels::removeGeometry(std::string const& geo_name,
                               GeoLib::GEOTYPE const type)
{
    if (type == GeoLib::GEOTYPE::SURFACE)
    {
        geo_objects_.removeSurfaceVec(geo_name);
    }
    if (type == GeoLib::GEOTYPE::POLYLINE)
    {
        geo_objects_.removePolylineVec(geo_name);
    }
    if (type == GeoLib::GEOTYPE::POINT)
    {
        geo_objects_.removePointVec(geo_name);
    }
}

void GEOModels::addPointVec(std::string const& name)
{
    geoModel_->addPointList(QString::fromStdString(name),
                            *geo_objects_.getPointVecObj(name));
    emit geoDataAdded(geoModel_, name, GeoLib::GEOTYPE::POINT);
}

void GEOModels::removePointVec(std::string const& name)
{
    assert(!geo_objects_.isPntVecUsed(name));

    emit geoDataRemoved(geoModel_, name, GeoLib::GEOTYPE::POINT);
    this->geoModel_->removeGeoList(name, GeoLib::GEOTYPE::POINT);
}

void GEOModels::addStationVec(std::string const& name)
{
    stationModel_->addStationList(QString::fromStdString(name),
                                  geo_objects_.getStationVec(name));
    emit stationVectorAdded(stationModel_, name);
}

void GEOModels::removeStationVec(std::string const& name)
{
    emit stationVectorRemoved(stationModel_, name);
    stationModel_->removeStationList(name);
}

void GEOModels::addPolylineVec(std::string const& name)
{
    geoModel_->addPolylineList(QString::fromStdString(name),
                               *geo_objects_.getPolylineVecObj(name));
    emit geoDataAdded(geoModel_, name, GeoLib::GEOTYPE::POLYLINE);
}

void GEOModels::appendPolylineVec(std::string const& name)
{
    this->geoModel_->appendPolylines(name,
                                     *geo_objects_.getPolylineVecObj(name));
}

void GEOModels::removePolylineVec(std::string const& name)
{
    emit geoDataRemoved(geoModel_, name, GeoLib::GEOTYPE::POLYLINE);
    this->geoModel_->removeGeoList(name, GeoLib::GEOTYPE::POLYLINE);
}

void GEOModels::addSurfaceVec(std::string const& name)
{
    geoModel_->addSurfaceList(QString::fromStdString(name),
                              *geo_objects_.getSurfaceVecObj(name));
    emit geoDataAdded(geoModel_, name, GeoLib::GEOTYPE::SURFACE);
}

void GEOModels::appendSurfaceVec(std::string const& name)
{
    geoModel_->appendSurfaces(name, *geo_objects_.getSurfaceVecObj(name));
}

void GEOModels::removeSurfaceVec(std::string const& name)
{
    emit geoDataRemoved(geoModel_, name, GeoLib::GEOTYPE::SURFACE);
    geoModel_->removeGeoList(name, GeoLib::GEOTYPE::SURFACE);
}

void GEOModels::renameGeometry(std::string const& old_name,
                               std::string const& new_name)
{
    geoModel_->renameGeometry(old_name, new_name);
    updateGeometry(new_name);
}

void GEOModels::connectPolylineSegments(
    const std::string& geoName, std::vector<std::size_t> const& indexlist,
    double const proximity, std::string const& ply_name, bool const closePly,
    bool const triangulatePly, std::string const& gmsh_path)
{
    GeoLib::PolylineVec* plyVec = geo_objects_.getPolylineVecObj(geoName);

    if (plyVec)
    {
        std::vector<GeoLib::Polyline*> const& polylines = *plyVec->getVector();
        std::vector<GeoLib::Polyline*> ply_list;
        std::transform(indexlist.begin(), indexlist.end(),
                       std::back_inserter(ply_list),
                       [polylines](auto const& ply_index) {
                           return polylines[ply_index];
                       });

        // connect polylines
        GeoLib::Polyline* new_line =
            GeoLib::Polyline::constructPolylineFromSegments(ply_list,
                                                            proximity);

        if (new_line)
        {
            // insert result in a new vector of polylines (because the GEOObjects::appendPolylines()-method requires a vector)
            std::vector<GeoLib::Polyline*> connected_ply;

            connected_ply.push_back(new_line);
            geo_objects_.appendPolylineVec(connected_ply, geoName);

            if (closePly)
            {
                new_line->closePolyline();

                if (triangulatePly)
                {
                    INFO(
                        "Creating a surface by triangulation of the polyline "
                        "...");
                    if (FileIO::createSurface(*new_line, geo_objects_, geoName,
                                              gmsh_path))
                    {
                        INFO("\t done");
                    }
                    else
                    {
                        WARN(
                            "\t Creating a surface by triangulation of the "
                            "polyline failed.");
                    }
                    plyVec = geo_objects_.getPolylineVecObj(geoName);
                }
            }

            if (!ply_name.empty())
            {
                plyVec->setNameOfElementByID(polylines.size(), ply_name);
            }
        }
        else
        {
            OGSError::box("Error connecting polyines.");
        }
    }
    else
    {
        OGSError::box("Corresponding geometry not found.");
    }
}

void GEOModels::addNameForElement(std::string const& geometry_name,
                                  GeoLib::GEOTYPE const object_type,
                                  std::size_t const id,
                                  std::string const& new_name)
{
    if (object_type == GeoLib::GEOTYPE::POINT)
    {
        geo_objects_.getPointVecObj(geometry_name)->setNameForElement(id, new_name);
    }
    else if (object_type == GeoLib::GEOTYPE::POLYLINE)
    {
        geo_objects_.getPolylineVecObj(geometry_name)->setNameForElement(id, new_name);
    }
    else if (object_type == GeoLib::GEOTYPE::SURFACE)
    {
        geo_objects_.getSurfaceVecObj(geometry_name)->setNameForElement(id, new_name);
    }
    else
        ERR("GEOModels::addNameForElement() - Unknown GEOTYPE {:s}.",
            GeoLib::convertGeoTypeToString(object_type));
}

void GEOModels::addNameForObjectPoints(const std::string &geometry_name,
                                       const GeoLib::GEOTYPE object_type,
                                       const std::string &geo_object_name,
                                       const std::string &new_name)
{
    const GeoLib::GeoObject* obj = geo_objects_.getGeoObject(geometry_name,
                                                      object_type,
                                                      geo_object_name);
    GeoLib::PointVec* pnt_vec = geo_objects_.getPointVecObj(geometry_name);
    if (object_type == GeoLib::GEOTYPE::POLYLINE)
    {
        const auto* ply = dynamic_cast<const GeoLib::Polyline*>(obj);
        std::size_t nPoints = ply->getNumberOfPoints();
        for (std::size_t i = 0; i < nPoints; i++)
        {
            pnt_vec->setNameForElement(
                ply->getPointID(i),
                new_name + "_Point" + std::to_string(ply->getPointID(i)));
        }
    }
    else if (object_type == GeoLib::GEOTYPE::SURFACE)
    {
        const auto* sfc = dynamic_cast<const GeoLib::Surface*>(obj);
        std::size_t nTriangles = sfc->getNumberOfTriangles();
        for (std::size_t i = 0; i < nTriangles; i++)
        {
            const GeoLib::Triangle* tri = (*sfc)[i];
            pnt_vec->setNameForElement(
                (*tri)[0], new_name + "_Point" + std::to_string((*tri)[0]));
            pnt_vec->setNameForElement(
                (*tri)[1], new_name + "_Point" + std::to_string((*tri)[1]));
            pnt_vec->setNameForElement(
                (*tri)[2], new_name + "_Point" + std::to_string((*tri)[2]));
        }
    }
    else
        ERR("GEOModels::addNameForObjectPoints() - Unknown GEOTYPE {:s}.",
            GeoLib::convertGeoTypeToString(object_type));
}
