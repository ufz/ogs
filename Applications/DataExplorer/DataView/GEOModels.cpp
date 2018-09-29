/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-09
 * \brief  Implementation of the GEOModels class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "GEOModels.h"

#include <logog/include/logog.hpp>

#include "Applications/FileIO/Legacy/createSurface.h"
#include "GeoTreeModel.h"
#include "StationTreeModel.h"

#include "OGSError.h"

#include "GeoLib/Triangle.h"

GEOModels::GEOModels(GeoLib::GEOObjects& geo_objects, QObject* parent /*= 0*/)
    : QObject(parent), _geo_objects(geo_objects)
{
    _geo_objects._callbacks = std::make_unique<GEOModelsCallbacks>(*this);

    _geoModel = new GeoTreeModel();
    _stationModel = new StationTreeModel();
}

GEOModels::~GEOModels()
{
    delete _stationModel;
    delete _geoModel;
}

void GEOModels::updateGeometry(const std::string &geo_name)
{
    GeoLib::PointVec* points (_geo_objects.getPointVecObj(geo_name));
    GeoLib::PolylineVec* lines (_geo_objects.getPolylineVecObj(geo_name));
    GeoLib::SurfaceVec* surfaces (_geo_objects.getSurfaceVecObj(geo_name));

    if (points)
    {
        emit geoDataRemoved(_geoModel, geo_name, GeoLib::GEOTYPE::POINT);
        this->_geoModel->removeGeoList(geo_name, GeoLib::GEOTYPE::POINT);
        _geoModel->addPointList(QString::fromStdString(geo_name), *points);
        emit geoDataAdded(_geoModel, geo_name, GeoLib::GEOTYPE::POINT);

        if (lines)
        {
            emit geoDataRemoved(_geoModel, geo_name, GeoLib::GEOTYPE::POLYLINE);
            this->_geoModel->removeGeoList(geo_name, GeoLib::GEOTYPE::POLYLINE);
            _geoModel->addPolylineList(QString::fromStdString(geo_name), *lines);
            emit geoDataAdded(_geoModel, geo_name, GeoLib::GEOTYPE::POLYLINE);
        }

        if (surfaces)
        {
            emit geoDataRemoved(_geoModel, geo_name, GeoLib::GEOTYPE::SURFACE);
            this->_geoModel->removeGeoList(geo_name, GeoLib::GEOTYPE::SURFACE);
            _geoModel->addSurfaceList(QString::fromStdString(geo_name), *surfaces);
            emit geoDataAdded(_geoModel, geo_name, GeoLib::GEOTYPE::SURFACE);
        }
    }
    else
        ERR("GEOModels::updateGeometry() - Geometry \"%s\" not found.", geo_name.c_str());
}

void GEOModels::removeGeometry(std::string const& geo_name,
                               GeoLib::GEOTYPE const type)
{
    if (type == GeoLib::GEOTYPE::SURFACE)
        _geo_objects.removeSurfaceVec(geo_name);
    if (type == GeoLib::GEOTYPE::POLYLINE)
        _geo_objects.removePolylineVec(geo_name);
    if (type == GeoLib::GEOTYPE::POINT)
        _geo_objects.removePointVec(geo_name);
}

void GEOModels::addPointVec(std::string const& name)
{
    _geoModel->addPointList(QString::fromStdString(name),
                            *_geo_objects.getPointVecObj(name));
    emit geoDataAdded(_geoModel, name, GeoLib::GEOTYPE::POINT);
}

void GEOModels::removePointVec(std::string const& name)
{
    assert(!_geo_objects.isPntVecUsed(name));

    emit geoDataRemoved(_geoModel, name, GeoLib::GEOTYPE::POINT);
    this->_geoModel->removeGeoList(name, GeoLib::GEOTYPE::POINT);
}

void GEOModels::addStationVec(std::string const& name)
{
    _stationModel->addStationList(QString::fromStdString(name),
                                  _geo_objects.getStationVec(name));
    emit stationVectorAdded(_stationModel, name);
}

void GEOModels::removeStationVec(std::string const& name)
{
    emit stationVectorRemoved(_stationModel, name);
    _stationModel->removeStationList(name);
}

void GEOModels::addPolylineVec(std::string const& name)
{
    _geoModel->addPolylineList(QString::fromStdString(name),
                               *_geo_objects.getPolylineVecObj(name));
    emit geoDataAdded(_geoModel, name, GeoLib::GEOTYPE::POLYLINE);
}

void GEOModels::appendPolylineVec(std::string const& name)
{
    this->_geoModel->appendPolylines(name,
                                     *_geo_objects.getPolylineVecObj(name));
}

void GEOModels::removePolylineVec(std::string const& name)
{
    emit geoDataRemoved(_geoModel, name, GeoLib::GEOTYPE::POLYLINE);
    this->_geoModel->removeGeoList(name, GeoLib::GEOTYPE::POLYLINE);
}

void GEOModels::addSurfaceVec(std::string const& name)
{
    _geoModel->addSurfaceList(QString::fromStdString(name),
                              *_geo_objects.getSurfaceVecObj(name));
    emit geoDataAdded(_geoModel, name, GeoLib::GEOTYPE::SURFACE);
}

void GEOModels::appendSurfaceVec(std::string const& name)
{
    _geoModel->appendSurfaces(name, *_geo_objects.getSurfaceVecObj(name));
}

void GEOModels::removeSurfaceVec(std::string const& name)
{
    emit geoDataRemoved(_geoModel, name, GeoLib::GEOTYPE::SURFACE);
    _geoModel->removeGeoList(name, GeoLib::GEOTYPE::SURFACE);
}

void GEOModels::renameGeometry(std::string const& old_name,
                               std::string const& new_name)
{
    _geoModel->renameGeometry(old_name, new_name);
    updateGeometry(new_name);
}

void GEOModels::connectPolylineSegments(
    const std::string& geoName,
    std::vector<std::size_t> const& indexlist,
    double const proximity,
    std::string const& ply_name,
    bool const closePly,
    bool const triangulatePly)
{
    GeoLib::PolylineVec* plyVec = _geo_objects.getPolylineVecObj(geoName);

    if (plyVec)
    {
        const std::vector<GeoLib::Polyline*>* polylines = plyVec->getVector();
        std::vector<GeoLib::Polyline*> ply_list;
        for (auto & elem : indexlist)
            ply_list.push_back( (*polylines)[elem] );

        // connect polylines
        GeoLib::Polyline* new_line = GeoLib::Polyline::constructPolylineFromSegments(
                ply_list,
                proximity);

        if (new_line)
        {
            // insert result in a new vector of polylines (because the GEOObjects::appendPolylines()-method requires a vector)
            std::vector<GeoLib::Polyline*> connected_ply;

            connected_ply.push_back(new_line);
            _geo_objects.appendPolylineVec(connected_ply, geoName);

            if (closePly)
            {
                new_line->closePolyline();

                if (triangulatePly)
                {
                    INFO(
                        "Creating a surface by triangulation of the polyline "
                        "...");
                    if (FileIO::createSurface(*new_line, _geo_objects, geoName))
                    {
                        INFO("\t done");
                    }
                    else
                    {
                        WARN(
                            "\t Creating a surface by triangulation of the "
                            "polyline failed.");
                    }
                    plyVec = _geo_objects.getPolylineVecObj(geoName);
                }
            }

            if (!ply_name.empty())
                plyVec->setNameOfElementByID(polylines->size(), ply_name);
        }
        else
            OGSError::box("Error connecting polyines.");
    }
    else
        OGSError::box("Corresponding geometry not found.");
}

void GEOModels::addNameForElement(std::string const& geometry_name,
                                  GeoLib::GEOTYPE const object_type,
                                  std::size_t const id,
                                  std::string const& new_name)
{
    if (object_type == GeoLib::GEOTYPE::POINT)
        _geo_objects.getPointVecObj(geometry_name)->setNameForElement(id, new_name);
    else if (object_type == GeoLib::GEOTYPE::POLYLINE)
        _geo_objects.getPolylineVecObj(geometry_name)->setNameForElement(id, new_name);
    else if (object_type == GeoLib::GEOTYPE::SURFACE)
        _geo_objects.getSurfaceVecObj(geometry_name)->setNameForElement(id, new_name);
    else
        ERR("GEOModels::addNameForElement() - Unknown GEOTYPE %s.",
            GeoLib::convertGeoTypeToString(object_type).c_str());
}

void GEOModels::addNameForObjectPoints(const std::string &geometry_name,
                                       const GeoLib::GEOTYPE object_type,
                                       const std::string &geo_object_name,
                                       const std::string &new_name)
{
    const GeoLib::GeoObject* obj = _geo_objects.getGeoObject(geometry_name,
                                                      object_type,
                                                      geo_object_name);
    GeoLib::PointVec* pnt_vec = _geo_objects.getPointVecObj(geometry_name);
    if (object_type == GeoLib::GEOTYPE::POLYLINE)
    {
        const auto* ply = dynamic_cast<const GeoLib::Polyline*>(obj);
        std::size_t nPoints = ply->getNumberOfPoints();
        for (std::size_t i = 0; i < nPoints; i++)
            pnt_vec->setNameForElement(
                ply->getPointID(i),
                new_name + "_Point" + std::to_string(ply->getPointID(i)));
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
        ERR("GEOModels::addNameForObjectPoints() - Unknown GEOTYPE %s.",
            GeoLib::convertGeoTypeToString(object_type).c_str());
}
