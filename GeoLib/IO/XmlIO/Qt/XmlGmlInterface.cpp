/**
 * \file
 * \author Karsten Rink
 * \date   2011-11-23
 * \brief  Implementation of the XmlGmlInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "XmlGmlInterface.h"

#include <QFile>
#include <QTextCodec>
#include <QtXml/QDomDocument>

#include <logog/include/logog.hpp>

#include "BaseLib/BuildInfo.h"
#include "BaseLib/FileFinder.h"
#include "GeoLib/Triangle.h"

namespace
{
void deletePolylines(std::unique_ptr<std::vector<GeoLib::Polyline*>> polylines)
{
    for (GeoLib::Polyline* line : *polylines)
        delete line;
}

void deleteSurfaces(std::unique_ptr<std::vector<GeoLib::Surface*>> surfaces)
{
    for (GeoLib::Surface* surface : *surfaces)
        delete surface;
}
void deleteGeometry(std::unique_ptr<std::vector<GeoLib::Point*>> points,
                    std::unique_ptr<std::vector<GeoLib::Polyline*>>
                        polylines,
                    std::unique_ptr<std::vector<GeoLib::Surface*>>
                        surfaces)
{
    for (GeoLib::Point* point : *points)
        delete point;
    deletePolylines(std::move(polylines));
    deleteSurfaces(std::move(surfaces));
}
}

namespace GeoLib
{
namespace IO
{
XmlGmlInterface::XmlGmlInterface(GeoLib::GEOObjects& geo_objs) :
XMLInterface(), XMLQtInterface("OpenGeoSysGLI.xsd"), _geo_objs(geo_objs)
{
}

int XmlGmlInterface::readFile(const QString &fileName)
{
    if(XMLQtInterface::readFile(fileName) == 0)
        return 0;

    QDomDocument doc("OGS-GLI-DOM");
    doc.setContent(_fileData);
    QDomElement docElement = doc.documentElement(); //OpenGeoSysGLI
    if (docElement.nodeName().compare("OpenGeoSysGLI"))
    {
        ERR("XmlGmlInterface::readFile() - Unexpected XML root.");
        return 0;
    }

    std::string gliName("[NN]");

    auto points = std::make_unique<std::vector<GeoLib::Point*>>();
    auto polylines = std::make_unique<std::vector<GeoLib::Polyline*>>();
    auto surfaces = std::make_unique<std::vector<GeoLib::Surface*>>();

    using MapNameId = std::map<std::string, std::size_t>;
    auto pnt_names = std::make_unique<MapNameId>();
    auto ply_names = std::make_unique<MapNameId>();
    auto sfc_names = std::make_unique<MapNameId>();

    QDomNodeList geoTypes = docElement.childNodes();
    for (int i = 0; i < geoTypes.count(); i++)
    {
        const QDomNode type_node(geoTypes.at(i));
        const QString nodeName = type_node.nodeName();
        if (nodeName.compare("name") == 0)
            if (type_node.toElement().text().isEmpty())
            {
                ERR("XmlGmlInterface::readFile(): <name>-tag is empty.")
                deleteGeometry(std::move(points), std::move(polylines),
                               std::move(surfaces));
                return 0;
            }
            else
                gliName = type_node.toElement().text().toStdString();
        else if (nodeName.compare("points") == 0)
        {
            readPoints(type_node, points.get(), pnt_names.get());

            // if names-map is empty, set it to nullptr because it is not needed
            if (pnt_names->empty())
            {
                pnt_names.reset(nullptr);
            }
            _geo_objs.addPointVec(std::move(points), gliName,
                                  std::move(pnt_names));
        }
        else if (nodeName.compare("polylines") == 0)
        {
            try
            {
                readPolylines(type_node, polylines.get(),
                              *_geo_objs.getPointVec(gliName),
                              _geo_objs.getPointVecObj(gliName)->getIDMap(),
                              ply_names.get());
            }
            catch (std::runtime_error const& err)
            {
                // further reading is aborted and it is necessary to clean up
                _geo_objs.removePointVec(gliName);
                throw;
            }
            // if names-map is empty, set it to nullptr because it is not needed
            if (ply_names->empty())
            {
                ply_names.reset(nullptr);
            }
        }
        else if (nodeName.compare("surfaces") == 0)
        {
            try
            {
                readSurfaces(type_node, surfaces.get(),
                             *_geo_objs.getPointVec(gliName),
                             _geo_objs.getPointVecObj(gliName)->getIDMap(),
                             sfc_names.get());
            }
            catch (std::runtime_error const& err)
            {
                // further reading is aborted and it is necessary to clean up
                _geo_objs.removePointVec(gliName);
                _geo_objs.removePolylineVec(gliName);
                throw;
            }

            // if names-map is empty, set it to nullptr because it is not needed
            if (sfc_names->empty())
            {
                sfc_names.reset(nullptr);
            }
        }
    }

    if (!polylines->empty())
        _geo_objs.addPolylineVec(std::move(polylines), gliName,
                                 std::move(ply_names));

    if (!surfaces->empty())
        _geo_objs.addSurfaceVec(std::move(surfaces), gliName,
                                std::move(sfc_names));
    return 1;
}

void XmlGmlInterface::readPoints(const QDomNode& pointsRoot,
                                 std::vector<GeoLib::Point*>* points,
                                 std::map<std::string, std::size_t>* pnt_names)
{
    char* pEnd;
    QDomElement point = pointsRoot.firstChildElement();
    while (!point.isNull())
    {
        _idx_map.insert (std::pair<std::size_t, std::size_t>(
            strtol((point.attribute("id")).toStdString().c_str(), &pEnd, 10), points->size()));
        GeoLib::Point* p = new GeoLib::Point(point.attribute("x").toDouble(),
                                             point.attribute("y").toDouble(),
                                             point.attribute("z").toDouble(),
                                             point.attribute("id").toInt());
        if (point.hasAttribute("name"))
            pnt_names->insert( std::pair<std::string, std::size_t>(
                                        point.attribute("name").toStdString(), points->size()) );

        points->push_back(p);
        point = point.nextSiblingElement();
    }
}

void XmlGmlInterface::readPolylines(
    const QDomNode& polylinesRoot,
    std::vector<GeoLib::Polyline*>* polylines,
    std::vector<GeoLib::Point*> const& points,
    const std::vector<std::size_t>& pnt_id_map,
    std::map<std::string, std::size_t>* ply_names)
{
    std::size_t idx(0);
    QDomElement polyline = polylinesRoot.firstChildElement();
    while (!polyline.isNull())
    {
        idx = polylines->size();
        polylines->push_back(new GeoLib::Polyline(points));

        if (polyline.hasAttribute("name")) {
            std::string const ply_name(
                polyline.attribute("name").toStdString()
            );
            std::map<std::string, std::size_t>::const_iterator it(
                ply_names->find(ply_name)
            );
            if (it == ply_names->end()) {
                ply_names->insert(std::pair<std::string, std::size_t>(ply_name, idx));
            } else {
                WARN("Polyline \"%s\" exists already. The polyline will be "
                    "inserted without a name.", ply_name.c_str());
            }
        }

        QDomElement point = polyline.firstChildElement();
        auto accessOrError =
            [this, &polyline](auto pt_idx) {
                auto search = _idx_map.find(pt_idx);
                if (search == _idx_map.end())
                {
                    std::string polyline_name;
                    if (polyline.hasAttribute("name"))
                        polyline_name =
                            polyline.attribute("name").toStdString();
                    OGS_FATAL(
                        "Polyline `%s' contains the point id `%d', but the "
                        "id is not in the point list.",
                        polyline_name.c_str(), pt_idx);
                }
                return search->second;
            };

        while (!point.isNull())
        {
            (*polylines)[idx]->addPoint(
                pnt_id_map[accessOrError(point.text().toInt())]);
            point = point.nextSiblingElement();
        }

        polyline = polyline.nextSiblingElement();
    }
}

void XmlGmlInterface::readSurfaces(
    const QDomNode& surfacesRoot,
    std::vector<GeoLib::Surface*>* surfaces,
    std::vector<GeoLib::Point*> const& points,
    const std::vector<std::size_t>& pnt_id_map,
    std::map<std::string, std::size_t>* sfc_names)
{
    QDomElement surface = surfacesRoot.firstChildElement();
    while (!surface.isNull())
    {
        surfaces->push_back(new GeoLib::Surface(points));

        if (surface.hasAttribute("name"))
            sfc_names->insert( std::pair<std::string, std::size_t>( surface.attribute("name").toStdString(),
                                                                    surfaces->size()-1) );

        auto accessOrError =
            [this, &surface](auto pt_idx) {
                auto search = _idx_map.find(pt_idx);
                if (search == _idx_map.end())
                {
                    std::string surface_name;
                    if (surface.hasAttribute("name"))
                        surface_name =
                            surface.attribute("name").toStdString();
                    OGS_FATAL(
                        "Surface `%s' contains the point id `%d', but the "
                        "id is not in the point list.",
                        surface_name.c_str(), pt_idx);
                }
                return search->second;
            };

        QDomElement element = surface.firstChildElement();
        while (!element.isNull())
        {
            std::size_t p1 =
                pnt_id_map[accessOrError(element.attribute("p1").toInt())];
            std::size_t p2 =
                pnt_id_map[accessOrError(element.attribute("p2").toInt())];
            std::size_t p3 =
                pnt_id_map[accessOrError(element.attribute("p3").toInt())];
            surfaces->back()->addTriangle(p1,p2,p3);
            element = element.nextSiblingElement();
        }

        surface = surface.nextSiblingElement();
    }
}
bool XmlGmlInterface::write()
{
    if (this->_exportName.empty())
    {
        ERR("XmlGmlInterface::write(): No geometry specified.");
        return false;
    }

    std::size_t nPoints = 0, nPolylines = 0, nSurfaces = 0;

    _out << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"; // xml definition
    //_out << "<?xml-stylesheet type=\"text/xsl\" href=\"OpenGeoSysGLI.xsl\"?>\n\n"; // stylefile definition

    QDomDocument doc("OGS-GML-DOM");
    QDomElement root = doc.createElement("OpenGeoSysGLI");
    root.setAttribute( "xmlns:ogs", "http://www.opengeosys.org" );
    root.setAttribute( "xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance" );

    doc.appendChild(root);

    QDomElement geoNameTag = doc.createElement("name");
    root.appendChild(geoNameTag);
    QDomText geoNameText = doc.createTextNode(QString::fromStdString(_exportName));
    geoNameTag.appendChild(geoNameText);

    // POINTS
    QDomElement pointsListTag = doc.createElement("points");
    root.appendChild(pointsListTag);

    const GeoLib::PointVec* pnt_vec (_geo_objs.getPointVecObj(_exportName));
    if (pnt_vec)
    {
        const std::vector<GeoLib::Point*>* points (pnt_vec->getVector());

        if (!points->empty())
        {
            nPoints = points->size();
            for (std::size_t i = 0; i < nPoints; i++)
            {
                QDomElement pointTag = doc.createElement("point");
                pointTag.setAttribute("id", QString::number(i));
                pointTag.setAttribute(
                    "x",
                    QString::number((*(*points)[i])[0], 'f',
                                    std::numeric_limits<double>::digits10));
                pointTag.setAttribute(
                    "y",
                    QString::number((*(*points)[i])[1], 'f',
                                    std::numeric_limits<double>::digits10));
                pointTag.setAttribute(
                    "z",
                    QString::number((*(*points)[i])[2], 'f',
                                    std::numeric_limits<double>::digits10));

                std::string const& point_name(pnt_vec->getItemNameByID(i));
                if (!point_name.empty())
                    pointTag.setAttribute("name",
                                          QString::fromStdString(point_name));

                pointsListTag.appendChild(pointTag);
            }
        }
        else
        {
            ERR("XmlGmlInterface::write(): Point vector empty, abort writing geometry.");
            return 0;
        }
    }
    else
    {
        ERR("XmlGmlInterface::write(): Did not found any point vector, abort writing geometry.");
        return 0;
    }

    // POLYLINES
    const GeoLib::PolylineVec* ply_vec (_geo_objs.getPolylineVecObj(_exportName));
    if (ply_vec)
    {
        const std::vector<GeoLib::Polyline*>* polylines (ply_vec->getVector());

        if (polylines)
        {
            if (!polylines->empty())
            {
                QDomElement plyListTag = doc.createElement("polylines");
                root.appendChild(plyListTag);
                nPolylines = polylines->size();
                for (std::size_t i = 0; i < nPolylines; i++)
                {
                    QDomElement polylineTag = doc.createElement("polyline");
                    polylineTag.setAttribute("id", QString::number(i));

                    std::string ply_name("");
                    if (ply_vec->getNameOfElementByID(i, ply_name)) {
                        polylineTag.setAttribute("name", QString::fromStdString(ply_name));
                    } else {
                        ply_name = std::to_string(i);
                        polylineTag.setAttribute("name", QString::fromStdString(ply_name));
                    }

                    plyListTag.appendChild(polylineTag);

                    nPoints = (*polylines)[i]->getNumberOfPoints();
                    for (std::size_t j = 0; j < nPoints; j++)
                    {
                        QDomElement plyPointTag = doc.createElement("pnt");
                        polylineTag.appendChild(plyPointTag);
                        QDomText plyPointText = doc.createTextNode(QString::number(((*polylines)[i])->getPointID(j)));
                        plyPointTag.appendChild(plyPointText);
                    }
                }
            }
            else
                INFO("XmlGmlInterface::write(): Polyline vector empty, no polylines written to file.");
        }
    }
    else
        INFO("XmlGmlInterface::write(): Did not found any polyline vector, no polylines written to file.");


    // SURFACES
    const GeoLib::SurfaceVec* sfc_vec (_geo_objs.getSurfaceVecObj(_exportName));
    if (sfc_vec)
    {
        const std::vector<GeoLib::Surface*>* surfaces (sfc_vec->getVector());

        if (surfaces)
        {
            if (!surfaces->empty())
            {
                QDomElement sfcListTag = doc.createElement("surfaces");
                root.appendChild(sfcListTag);
                nSurfaces = surfaces->size();
                for (std::size_t i = 0; i < nSurfaces; i++)
                {
                    QDomElement surfaceTag = doc.createElement("surface");
                    surfaceTag.setAttribute("id", QString::number(i));

                    std::string sfc_name("");
                    if (sfc_vec->getNameOfElementByID(i, sfc_name))
                        surfaceTag.setAttribute("name",
                                                QString::fromStdString(
                                                        sfc_name));

                    sfcListTag.appendChild(surfaceTag);

                    // writing the elements compromising the surface
                    std::size_t nElements = ((*surfaces)[i])->getNumberOfTriangles();
                    for (std::size_t j = 0; j < nElements; j++)
                    {
                        QDomElement elementTag = doc.createElement("element");
                        elementTag.setAttribute("p1", QString::number((*(*(*surfaces)[i])[j])[0]));
                        elementTag.setAttribute("p2", QString::number((*(*(*surfaces)[i])[j])[1]));
                        elementTag.setAttribute("p3", QString::number((*(*(*surfaces)[i])[j])[2]));
                        surfaceTag.appendChild(elementTag);
                    }
                }
            }
            else
                INFO("XmlGmlInterface::write(): Surface vector empty, no surfaces written to file.");
        }
    }
    else
        INFO("XmlGmlInterface::write(): Did not found any surface vector, no surfaces written to file.");


    //insertStyleFileDefinition(filename);
    std::string xml = doc.toString().toStdString();
    _out << xml;

    return true;
}

}
}
