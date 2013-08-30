/**
 * \file
 * \author Karsten Rink
 * \date   2010-01-25
 * \brief  Definition of the SHPInterface class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * @file SHPInterface.cpp
 * @date 25/01/2010
 * @author Karsten Rink
 */

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "SHPInterface.h"

// BaseLib
#include "StringTools.h"

// MathLib
#include "AnalyticalGeometry.h"

// MeshLib
#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"
#include "Elements/Tri.h"
#include "Elements/Quad.h"



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

void SHPInterface::readSHPFile(const std::string &filename, OGSType choice, std::string listName)
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
		std::vector<GeoLib::Point*>* points = new std::vector<GeoLib::Point*>();
		SHPObject* hSHPObject;

		for (int i = 0; i < numberOfElements; i++) {
			hSHPObject = SHPReadObject(hSHP, i);

			GeoLib::Point* pnt =
			        new GeoLib::Point(*(hSHPObject->padfX), *(hSHPObject->padfY),
			                          *(hSHPObject->padfZ));
			points->push_back(pnt);
		}

		_geoObjects->addPointVec(points, listName);
		SHPDestroyObject(hSHPObject); // de-allocate SHPObject
	}
}

void SHPInterface::readStations(const SHPHandle &hSHP, int numberOfElements, std::string listName)
{
	if (numberOfElements > 0) {
		std::vector<GeoLib::Point*>* stations(new std::vector<GeoLib::Point*>);
		stations->reserve(numberOfElements);
		SHPObject* hSHPObject;

		for (int i = 0; i < numberOfElements; i++) {
			hSHPObject = SHPReadObject(hSHP, i);
			GeoLib::Station* stn = GeoLib::Station::createStation(BaseLib::number2str(i),
			                                                      *(hSHPObject->padfX),
			                                                      *(hSHPObject->padfY),
			                                                      *(hSHPObject->padfZ));
			stations->push_back(stn);
		}

		_geoObjects->addStationVec(stations, listName);
		SHPDestroyObject(hSHPObject); // de-allocate SHPObject
	}
}

void SHPInterface::readPolylines(const SHPHandle &hSHP, int numberOfElements, std::string listName)
{
	std::size_t noOfPoints = 0, noOfParts = 0;
	std::vector<GeoLib::Point*>* points = new std::vector<GeoLib::Point*>();
	std::vector<GeoLib::Polyline*>* lines = new std::vector<GeoLib::Polyline*>();
	SHPObject* hSHPObject;

	// for each polyline)
	for (int i = 0; i < numberOfElements; i++) {
		hSHPObject = SHPReadObject(hSHP, i);
		noOfPoints = hSHPObject->nVertices;
		noOfParts = hSHPObject->nParts;

		for (std::size_t p = 0; p < noOfParts; p++) {
			int firstPnt = *(hSHPObject->panPartStart + p);
			int lastPnt = (p < (noOfParts - 1)) ? *(hSHPObject->panPartStart + p + 1) : noOfPoints;

			GeoLib::Polyline* line = new GeoLib::Polyline(*points);

			// for each point in that polyline
			for (int j = firstPnt; j < lastPnt; j++) {
				GeoLib::Point* pnt = new GeoLib::Point(*(hSHPObject->padfX + j),
				                                       *(hSHPObject->padfY + j),
				                                       *(hSHPObject->padfZ + j));
				points->push_back(pnt);
				line->addPoint(points->size() - 1);
			}

			// add polyline to polyline vector
			lines->push_back(line);
		}
	}

	if (numberOfElements > 0) {
		// add points vector to GEOObjects (and check for duplicate points)
		_geoObjects->addPointVec(points, listName);

		// adjust indeces of polylines, remove zero length elements and add vector to GEOObjects
		this->adjustPolylines(lines, _geoObjects->getPointVecObj(listName)->getIDMap());
		_geoObjects->addPolylineVec(lines, listName);

		SHPDestroyObject(hSHPObject); // de-allocate SHPObject
	}
}

void SHPInterface::readPolygons(const SHPHandle &hSHP, int numberOfElements, std::string listName)
{
	this->readPolylines(hSHP, numberOfElements, listName);

	const std::vector<GeoLib::Polyline*>* polylines(_geoObjects->getPolylineVec(listName));
	std::vector<GeoLib::Surface*>* sfc_vec(new std::vector<GeoLib::Surface*>);

	for (std::vector<GeoLib::Polyline*>::const_iterator poly_it(polylines->begin()); poly_it
	                != polylines->end(); poly_it++) {
		GeoLib::Surface* sfc(GeoLib::Surface::createSurface(*(*poly_it)));
		if (sfc)
			sfc_vec->push_back(sfc);
		else {
			WARN("SHPInterface::readPolygons(): Could not triangulate polygon.")
		}
	}

	if (!sfc_vec->empty())
		_geoObjects->addSurfaceVec(sfc_vec, listName);
}

void SHPInterface::adjustPolylines(std::vector<GeoLib::Polyline*>* lines,
                                   std::vector<std::size_t> id_map)

{
	for (std::size_t i = 0; i < lines->size(); i++) {
		GeoLib::Polyline* line((*lines)[i]);
		std::size_t previous_pnt_id(std::numeric_limits<std::size_t>::max());

		for (std::size_t j = 0; j < line->getNumberOfPoints(); j++) {
			std::size_t jth_pnt_id(id_map[line->getPointID(j)]);
			if (previous_pnt_id == jth_pnt_id) {
				line->removePoint(j);
				j--;
			} else {
				line->setPointID(j, jth_pnt_id);
			}
			previous_pnt_id = jth_pnt_id;
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

	unsigned nElements (mesh.getNElements());
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
		if ((e->getGeomType() == MeshElemType::TRIANGLE) || (e->getGeomType() == MeshElemType::QUAD))
		{
			// write element ID and material group to DBF-file
			DBFWriteIntegerAttribute(hDBF, polygon_id, elem_id_field, i);
			DBFWriteIntegerAttribute(hDBF, polygon_id, mat_field, e->getValue());

			unsigned nNodes (e->getNNodes());
			padfX = new double(nNodes+1);
			padfY = new double(nNodes+1);
			padfZ = new double(nNodes+1);
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

			SHPObject *object = SHPCreateObject(SHPT_POLYGON, polygon_id++, 0, 0, NULL, ++nNodes, padfX, padfY, padfZ, NULL);
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
