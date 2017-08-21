/**
 * \file
 * \author Karsten Rink
 * \date   2010-01-25
 * \brief  Definition of the SHPInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * @file SHPInterface.cpp
 * @date 25/01/2010
 * @author Karsten Rink
 */

#include "SHPInterface.h"

#include <logog/include/logog.hpp>

#include "Applications/FileIO/Legacy/createSurface.h"

#include "GeoLib/AnalyticalGeometry.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Polyline.h"
#include "GeoLib/Point.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Elements/Quad.h"

namespace FileIO
{

bool SHPInterface::readSHPInfo(const std::string &filename, int &shapeType, int &numberOfEntities)
{
    SHPHandle hSHP = SHPOpen(filename.c_str(), "rb");
    if (!hSHP)
        return false;

    double padfMinBound[4], padfMaxBound[4];

    // The SHPGetInfo() function retrieves various information about shapefile as a whole.
    // The bounds are read from the file header, and may be inaccurate if the file was improperly generated.
    SHPGetInfo(hSHP, &numberOfEntities, &shapeType, padfMinBound, padfMaxBound);

    SHPClose(hSHP);
    return true;
}

void SHPInterface::readSHPFile(const std::string &filename, OGSType choice, const std::string &listName)
{
    int shapeType, numberOfElements;
    double padfMinBound[4], padfMaxBound[4];

    SHPHandle hSHP = SHPOpen(filename.c_str(), "rb");
    SHPGetInfo(hSHP, &numberOfElements, &shapeType, padfMinBound, padfMaxBound);

    if (((shapeType - 1) % 10 == 0) && (choice == SHPInterface::OGSType::POINT))
        readPoints(hSHP, numberOfElements, listName);
    if (((shapeType - 1) % 10 == 0) && (choice == SHPInterface::OGSType::STATION))
        readStations(hSHP, numberOfElements, listName);
    if (((shapeType - 3) % 10 == 0 || (shapeType - 5) % 10 == 0) && (choice
                    == SHPInterface::OGSType::POLYLINE))
        readPolylines(hSHP, numberOfElements, listName);
    if (((shapeType - 3) % 10 == 0 || (shapeType - 5) % 10 == 0) && (choice
                    == SHPInterface::OGSType::POLYGON))
        readPolygons(hSHP, numberOfElements, listName);
}

void SHPInterface::readPoints(const SHPHandle &hSHP, int numberOfElements, std::string listName)
{
    if (numberOfElements > 0) {
        auto points = std::make_unique<std::vector<GeoLib::Point*>>();
        SHPObject* hSHPObject;

        for (int i = 0; i < numberOfElements; i++) {
            hSHPObject = SHPReadObject(hSHP, i);

            auto* pnt =
                new GeoLib::Point(*(hSHPObject->padfX), *(hSHPObject->padfY),
                                  *(hSHPObject->padfZ));
            points->push_back(pnt);
        }

        _geoObjects.addPointVec(std::move(points), listName);
        SHPDestroyObject(hSHPObject); // de-allocate SHPObject
    }
}

void SHPInterface::readStations(const SHPHandle &hSHP, int numberOfElements, std::string listName)
{
    if (numberOfElements > 0) {
        auto stations = std::make_unique<std::vector<GeoLib::Point*>>();
        stations->reserve(numberOfElements);
        SHPObject* hSHPObject;

        for (int i = 0; i < numberOfElements; i++) {
            hSHPObject = SHPReadObject(hSHP, i);
            GeoLib::Station* stn = GeoLib::Station::createStation(std::to_string(i),
                                                                  *(hSHPObject->padfX),
                                                                  *(hSHPObject->padfY),
                                                                  *(hSHPObject->padfZ));
            stations->push_back(stn);
        }

        _geoObjects.addStationVec(std::move(stations), listName);
        SHPDestroyObject(hSHPObject); // de-allocate SHPObject
    }
}

void SHPInterface::readPolylines(const SHPHandle &hSHP, int numberOfElements, std::string listName)
{
    if (numberOfElements <= 0)
        return;
    auto pnts = std::make_unique<std::vector<GeoLib::Point*>>();
    auto lines = std::make_unique<std::vector<GeoLib::Polyline*>>();

    std::size_t pnt_id(0);
    // for each polyline
    for (int i = 0; i < numberOfElements; ++i) {
        SHPObject *hSHPObject = SHPReadObject(hSHP, i);
        int const noOfPoints = hSHPObject->nVertices;
        int const noOfParts = hSHPObject->nParts;

        for (int p = 0; p < noOfParts; ++p) {
            int const firstPnt = *(hSHPObject->panPartStart + p);
            int const lastPnt = (p<(noOfParts - 1)) ?
                *(hSHPObject->panPartStart + p + 1) : noOfPoints;

            // for each point in that polyline
            for (int j = firstPnt; j < lastPnt; ++j) {
                pnts->push_back(new GeoLib::Point(*(hSHPObject->padfX+j),
                    *(hSHPObject->padfY+j), *(hSHPObject->padfZ+j), pnt_id));
                pnt_id++;
            }
        }
        SHPDestroyObject(hSHPObject); // de-allocate SHPObject
    }

    _geoObjects.addPointVec(std::move(pnts), listName);
    GeoLib::PointVec const& points(*(_geoObjects.getPointVecObj(listName)));
    std::vector<std::size_t> const& pnt_id_map(points.getIDMap());

    pnt_id = 0;
    for (int i = 0; i < numberOfElements; ++i) {
        SHPObject* hSHPObject = SHPReadObject(hSHP, i);
        int const noOfPoints = hSHPObject->nVertices;
        int const noOfParts = hSHPObject->nParts;

        for (int p = 0; p < noOfParts; ++p) {
            int const firstPnt = *(hSHPObject->panPartStart + p);
            int const lastPnt = (p<(noOfParts - 1)) ?
                *(hSHPObject->panPartStart + p + 1) : noOfPoints;

            auto* line = new GeoLib::Polyline(*points.getVector());

            // create polyline
            for (int j = firstPnt; j < lastPnt; ++j) {
                line->addPoint(pnt_id_map[pnt_id]);
                pnt_id++;
            }
            // add polyline to polyline vector
            lines->push_back(line);
        }
        SHPDestroyObject(hSHPObject); // de-allocate SHPObject
    }
    _geoObjects.addPolylineVec(std::move(lines), listName);
}

void SHPInterface::readPolygons(const SHPHandle &hSHP, int numberOfElements, const std::string &listName)
{
    readPolylines(hSHP, numberOfElements, listName);

    auto const polylines = _geoObjects.getPolylineVec(listName);
    auto sfc_vec = std::make_unique<std::vector<GeoLib::Surface*>>();

    for (auto const* polyline : *polylines)
    {
        INFO("Creating a surface by triangulation of the polyline ...");
        if (FileIO::createSurface(*polyline, _geoObjects, listName))
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

bool SHPInterface::write2dMeshToSHP(const std::string &file_name, const MeshLib::Mesh &mesh)
{
    if (mesh.getDimension()!=2)
    {
        ERR ("SHPInterface::write2dMeshToSHP(): Mesh to Shape conversion is only working for 2D Meshes.");
        return false;
    }

    unsigned nElements (mesh.getNumberOfElements());
    if (nElements<1)
    {
        ERR ("SHPInterface::write2dMeshToSHP(): Mesh contains no elements.");
        return false;
    }

    if (nElements>10E+7) // DBF-export requires a limit, 10 mio seems good for now
    {
        ERR ("SHPInterface::write2dMeshToSHP(): Mesh contains too many elements for currently implemented DBF-boundaries.");
        return false;
    }

    SHPHandle hSHP = SHPCreate(file_name.c_str(), SHPT_POLYGON);
    DBFHandle hDBF = DBFCreate(file_name.c_str());
    int elem_id_field = DBFAddField(hDBF, "Elem_ID", FTInteger, 7, 0); // allows integers of length "7", i.e. 10mio-1 elements
    int mat_field = DBFAddField(hDBF, "Material", FTInteger, 7, 0);
    int node0_field = DBFAddField(hDBF, "Node0", FTInteger, 7, 0);
    int node1_field = DBFAddField(hDBF, "Node1", FTInteger, 7, 0);
    int node2_field = DBFAddField(hDBF, "Node2", FTInteger, 7, 0);

    unsigned polygon_id (0);
    double* padfX;
    double* padfY;
    double* padfZ;
    for (unsigned i=0; i<nElements; ++i)
    {
        const MeshLib::Element* e (mesh.getElement(i));

        // ignore all elements except triangles and quads
        if ((e->getGeomType() == MeshLib::MeshElemType::TRIANGLE) ||
            (e->getGeomType() == MeshLib::MeshElemType::QUAD))
        {
            // write element ID and material group to DBF-file
            DBFWriteIntegerAttribute(hDBF, polygon_id, elem_id_field, i);
            if (mesh.getProperties().existsPropertyVector<int>("MaterialIDs"))
            {
                auto const* const materialIds =
                    mesh.getProperties().getPropertyVector<int>("MaterialIDs");
                DBFWriteIntegerAttribute(hDBF, polygon_id, mat_field, (*materialIds)[i]);
            }
            unsigned nNodes (e->getNumberOfBaseNodes());
            padfX = new double[nNodes+1];
            padfY = new double[nNodes+1];
            padfZ = new double[nNodes+1];
            for (unsigned j=0; j<nNodes; ++j)
            {
                padfX[j]=(*e->getNode(j))[0];
                padfY[j]=(*e->getNode(j))[1];
                padfZ[j]=(*e->getNode(j))[2];
            }
            // Last node == first node to close the polygon
            padfX[nNodes]=(*e->getNode(0))[0];
            padfY[nNodes]=(*e->getNode(0))[1];
            padfZ[nNodes]=(*e->getNode(0))[2];
            // write the first three node ids to the dbf-file (this also specifies a QUAD uniquely)
            DBFWriteIntegerAttribute(hDBF, polygon_id, node0_field, e->getNode(0)->getID());
            DBFWriteIntegerAttribute(hDBF, polygon_id, node1_field, e->getNode(1)->getID());
            DBFWriteIntegerAttribute(hDBF, polygon_id, node2_field, e->getNode(2)->getID());

            SHPObject* object =
                SHPCreateObject(SHPT_POLYGON, polygon_id++, 0, nullptr, nullptr,
                                ++nNodes, padfX, padfY, padfZ, nullptr);
            SHPWriteObject(hSHP, -1, object);

            // Note: cleaning up the coordinate arrays padfX, -Y, -Z results in a crash, I assume that shapelib removes them
            delete object;
        }
    }

    SHPClose(hSHP);
    DBFClose(hDBF);
    INFO ("Shape export of 2D mesh \"%s\" successful.", mesh.getName().c_str());

    return true;
}

}
