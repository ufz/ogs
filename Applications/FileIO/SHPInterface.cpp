/**
 * \file
 * \author Karsten Rink
 * \date   2010-01-25
 * \brief  Definition of the SHPInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * @file SHPInterface.cpp
 * @date 25/01/2010
 * @author Karsten Rink
 */

#include "SHPInterface.h"

#include "Applications/FileIO/Legacy/createSurface.h"
#include "BaseLib/Logging.h"
#include "GeoLib/AnalyticalGeometry.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Point.h"
#include "GeoLib/Polyline.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

namespace FileIO
{
bool SHPInterface::readSHPInfo(const std::string& filename, int& shapeType,
                               int& numberOfEntities)
{
    SHPHandle hSHP = SHPOpen(filename.c_str(), "rb");
    if (!hSHP)
    {
        return false;
    }

    double padfMinBound[4];
    double padfMaxBound[4];

    // The SHPGetInfo() function retrieves various information about shapefile
    // as a whole. The bounds are read from the file header, and may be
    // inaccurate if the file was improperly generated.
    SHPGetInfo(hSHP, &numberOfEntities, &shapeType, padfMinBound, padfMaxBound);

    SHPClose(hSHP);
    return true;
}

void SHPInterface::readSHPFile(const std::string& filename, OGSType choice,
                               const std::string& listName,
                               std::string const& gmsh_path)
{
    int shapeType;
    int numberOfElements;
    double padfMinBound[4];
    double padfMaxBound[4];

    SHPHandle hSHP = SHPOpen(filename.c_str(), "rb");
    SHPGetInfo(hSHP, &numberOfElements, &shapeType, padfMinBound, padfMaxBound);

    if (((shapeType - 1) % 10 == 0) && (choice == SHPInterface::OGSType::POINT))
    {
        readPoints(hSHP, numberOfElements, listName);
    }
    if (((shapeType - 1) % 10 == 0) &&
        (choice == SHPInterface::OGSType::STATION))
    {
        readStations(hSHP, numberOfElements, listName);
    }
    if (((shapeType - 3) % 10 == 0 || (shapeType - 5) % 10 == 0) &&
        (choice == SHPInterface::OGSType::POLYLINE))
    {
        readPolylines(hSHP, numberOfElements, listName);
    }
    if (((shapeType - 3) % 10 == 0 || (shapeType - 5) % 10 == 0) &&
        (choice == SHPInterface::OGSType::POLYGON))
    {
        readPolygons(hSHP, numberOfElements, listName, gmsh_path);
    }
}

void SHPInterface::readPoints(const SHPHandle& hSHP, int numberOfElements,
                              std::string listName)
{
    if (numberOfElements > 0)
    {
        std::vector<GeoLib::Point*> points;
        SHPObject* hSHPObject;

        for (int i = 0; i < numberOfElements; i++)
        {
            hSHPObject = SHPReadObject(hSHP, i);

            auto* pnt =
                new GeoLib::Point(*(hSHPObject->padfX), *(hSHPObject->padfY),
                                  *(hSHPObject->padfZ));
            points.push_back(pnt);
        }

        _geoObjects.addPointVec(std::move(points), listName,
                                GeoLib::PointVec::NameIdMap{});
        SHPDestroyObject(hSHPObject);  // de-allocate SHPObject
    }
}

void SHPInterface::readStations(const SHPHandle& hSHP, int numberOfElements,
                                std::string listName)
{
    if (numberOfElements > 0)
    {
        std::vector<GeoLib::Point*> stations;
        stations.reserve(numberOfElements);
        SHPObject* hSHPObject;

        for (int i = 0; i < numberOfElements; i++)
        {
            hSHPObject = SHPReadObject(hSHP, i);
            GeoLib::Station* stn =
                GeoLib::Station::createStation(std::to_string(i),
                                               *(hSHPObject->padfX),
                                               *(hSHPObject->padfY),
                                               *(hSHPObject->padfZ));
            stations.push_back(stn);
        }

        _geoObjects.addStationVec(std::move(stations), listName);
        SHPDestroyObject(hSHPObject);  // de-allocate SHPObject
    }
}

void SHPInterface::readPolylines(const SHPHandle& hSHP, int numberOfElements,
                                 std::string listName)
{
    if (numberOfElements <= 0)
    {
        return;
    }
    std::vector<GeoLib::Point*> pnts;
    std::vector<GeoLib::Polyline*> lines;

    std::size_t pnt_id(0);
    // for each polyline
    for (int i = 0; i < numberOfElements; ++i)
    {
        SHPObject* hSHPObject = SHPReadObject(hSHP, i);
        int const noOfPoints = hSHPObject->nVertices;
        int const noOfParts = hSHPObject->nParts;

        for (int p = 0; p < noOfParts; ++p)
        {
            int const firstPnt = *(hSHPObject->panPartStart + p);
            int const lastPnt = (p < (noOfParts - 1))
                                    ? *(hSHPObject->panPartStart + p + 1)
                                    : noOfPoints;

            // for each point in that polyline
            for (int j = firstPnt; j < lastPnt; ++j)
            {
                pnts.push_back(new GeoLib::Point(
                    *(hSHPObject->padfX + j), *(hSHPObject->padfY + j),
                    *(hSHPObject->padfZ + j), pnt_id));
                pnt_id++;
            }
        }
        SHPDestroyObject(hSHPObject);  // de-allocate SHPObject
    }

    _geoObjects.addPointVec(std::move(pnts), listName,
                            GeoLib::PointVec::NameIdMap{});
    GeoLib::PointVec const& points(*(_geoObjects.getPointVecObj(listName)));
    std::vector<std::size_t> const& pnt_id_map(points.getIDMap());

    pnt_id = 0;
    for (int i = 0; i < numberOfElements; ++i)
    {
        SHPObject* hSHPObject = SHPReadObject(hSHP, i);
        int const noOfPoints = hSHPObject->nVertices;
        int const noOfParts = hSHPObject->nParts;

        for (int p = 0; p < noOfParts; ++p)
        {
            // output for the first part of multipart polyline
            if (noOfParts > 1 && p == 0)
            {
                INFO(
                    "Polygon {:d} consists of {:d} parts (PolylineIDs "
                    "{:d}-{:d}).",
                    i, noOfParts, lines.size(), lines.size() + noOfParts - 1);
            }

            int const firstPnt = *(hSHPObject->panPartStart + p);
            int const lastPnt = (p < (noOfParts - 1))
                                    ? *(hSHPObject->panPartStart + p + 1)
                                    : noOfPoints;

            auto* line = new GeoLib::Polyline(*points.getVector());

            // create polyline
            for (int j = firstPnt; j < lastPnt; ++j)
            {
                line->addPoint(pnt_id_map[pnt_id]);
                pnt_id++;
            }
            // add polyline to polyline vector
            lines.push_back(line);
        }
        SHPDestroyObject(hSHPObject);  // de-allocate SHPObject
    }
    _geoObjects.addPolylineVec(std::move(lines), listName,
                               GeoLib::PolylineVec::NameIdMap{});
}

void SHPInterface::readPolygons(const SHPHandle& hSHP, int numberOfElements,
                                const std::string& listName,
                                std::string const& gmsh_path)
{
    readPolylines(hSHP, numberOfElements, listName);

    auto const polylines = _geoObjects.getPolylineVec(listName);

    for (auto const* polyline : *polylines)
    {
        INFO("Creating a surface by triangulation of the polyline ...");
        if (FileIO::createSurface(*polyline, _geoObjects, listName, gmsh_path))
        {
            INFO("\t done");
        }
        else
        {
            WARN(
                "\t Creating a surface by triangulation of the polyline "
                "failed.");
        }
    }
}

bool SHPInterface::write2dMeshToSHP(const std::string& file_name,
                                    const MeshLib::Mesh& mesh)
{
    if (mesh.getDimension() != 2)
    {
        ERR("SHPInterface::write2dMeshToSHP(): Mesh to Shape conversion is "
            "only working for 2D Meshes.");
        return false;
    }

    std::size_t const n_elements = mesh.getNumberOfElements();
    if (n_elements < 1)
    {
        ERR("SHPInterface::write2dMeshToSHP(): Mesh contains no elements.");
        return false;
    }

    // DBF-export requires a limit of records because of the limits to the
    // integer values. The exponent controls maximum number of elements and
    // the maximum number of digits written in the DBF file.
    std::size_t const max_exp = 8;
    if (n_elements >= std::pow(10, max_exp))
    {
        ERR("SHPInterface::write2dMeshToSHP(): Mesh contains too many elements "
            "for currently implemented DBF-boundaries.");
        return false;
    }

    SHPHandle hSHP = SHPCreate(file_name.c_str(), SHPT_POLYGON);
    DBFHandle hDBF = DBFCreate(file_name.c_str());
    int const elem_id_field =
        DBFAddField(hDBF, "Elem_ID", FTInteger, max_exp, 0);

    // Writing mesh elements to shape file
    std::size_t shp_record(0);
    std::vector<MeshLib::Element*> const& elems = mesh.getElements();
    for (MeshLib::Element const* const e : elems)
    {
        // ignore all elements except triangles and quads
        if ((e->getGeomType() == MeshLib::MeshElemType::TRIANGLE) ||
            (e->getGeomType() == MeshLib::MeshElemType::QUAD))
        {
            SHPObject* object = createShapeObject(*e, shp_record);
            SHPWriteObject(hSHP, -1, object);
            SHPDestroyObject(object);

            // write element ID to DBF-file
            DBFWriteIntegerAttribute(hDBF, shp_record, elem_id_field,
                                     e->getID());
            shp_record++;
        }
    }
    SHPClose(hSHP);

    // Write scalar arrays to database file
    int const n_recs = DBFGetRecordCount(hDBF);
    for (auto [name, property] : mesh.getProperties())
    {
        if (auto p = dynamic_cast<MeshLib::PropertyVector<int>*>(property))
        {
            int const field = DBFAddField(hDBF, name.c_str(), FTInteger, 16, 0);
            for (int i = 0; i < n_recs; ++i)
            {
                std::size_t const elem_idx =
                    DBFReadIntegerAttribute(hDBF, i, elem_id_field);
                DBFWriteIntegerAttribute(hDBF, i, field, (*p)[elem_idx]);
            }
        }
        else if (auto p =
                     dynamic_cast<MeshLib::PropertyVector<double>*>(property))
        {
            int const field = DBFAddField(hDBF, name.c_str(), FTDouble, 33, 16);
            for (int i = 0; i < n_recs; ++i)
            {
                std::size_t const elem_idx =
                    DBFReadIntegerAttribute(hDBF, i, elem_id_field);
                DBFWriteDoubleAttribute(hDBF, i, field, (*p)[elem_idx]);
            }
        }
    }
    DBFClose(hDBF);
    INFO("Shape export of 2D mesh '{:s}' finished.", mesh.getName());
    return true;
}

SHPObject* SHPInterface::createShapeObject(MeshLib::Element const& e,
                                           std::size_t const shp_record)
{
    unsigned const nNodes(e.getNumberOfBaseNodes());
    double* padfX = new double[nNodes + 1];
    double* padfY = new double[nNodes + 1];
    double* padfZ = new double[nNodes + 1];
    for (unsigned j = 0; j < nNodes; ++j)
    {
        padfX[j] = (*e.getNode(j))[0];
        padfY[j] = (*e.getNode(j))[1];
        padfZ[j] = (*e.getNode(j))[2];
    }
    // Last node == first node to close the polygon
    padfX[nNodes] = (*e.getNode(0))[0];
    padfY[nNodes] = (*e.getNode(0))[1];
    padfZ[nNodes] = (*e.getNode(0))[2];

    // the generated shape object now handles the memory for padfX/Y/Z
    return SHPCreateObject(SHPT_POLYGON, shp_record, 0, nullptr, nullptr,
                           nNodes + 1, padfX, padfY, padfZ, nullptr);
}

}  // namespace FileIO
