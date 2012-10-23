/**
 * \file SHPInterface.cpp
 * 25/01/2010 KR Initial implementation
 */

#include "SHPInterface.h"
#include "StringTools.h"

// MathLib
#include "AnalyticalGeometry.h"

bool SHPInterface::readSHPInfo(const std::string &filename, int &shapeType, int &numberOfEntities)
{
	SHPHandle hSHP = SHPOpen(filename.c_str(),"rb");
	if(!hSHP)
		return false;

	double padfMinBound[4], padfMaxBound[4];

	// The SHPGetInfo() function retrieves various information about shapefile as a whole.
	// The bounds are read from the file header, and may be inaccurate if the file was improperly generated.
	SHPGetInfo( hSHP, &numberOfEntities, &shapeType, padfMinBound, padfMaxBound );

	SHPClose(hSHP);
	return true;
}

void SHPInterface::readSHPFile(const std::string &filename, OGSType choice, std::string listName)
{
	int shapeType, numberOfElements;
	double padfMinBound[4], padfMaxBound[4];

	SHPHandle hSHP = SHPOpen(filename.c_str(),"rb");
	SHPGetInfo( hSHP, &numberOfElements, &shapeType, padfMinBound, padfMaxBound );

	if ( ((shapeType - 1) % 10 == 0)  &&  (choice == SHPInterface::POINT) )
		readPoints(hSHP, numberOfElements, listName);
	if ( ((shapeType - 1) % 10 == 0)  &&  (choice == SHPInterface::STATION) )
		readStations(hSHP, numberOfElements, listName);
	if ( ((shapeType - 3) % 10 == 0 ||
	      (shapeType - 5) % 10 == 0)  &&  (choice == SHPInterface::POLYLINE) )
		readPolylines(hSHP, numberOfElements, listName);
	if ( ((shapeType - 3) % 10 == 0 ||
	      (shapeType - 5) % 10 == 0)  &&  (choice == SHPInterface::POLYGON) )
		readPolygons(hSHP, numberOfElements, listName);
}

void SHPInterface::readPoints(const SHPHandle &hSHP, int numberOfElements, std::string listName)
{
	if (numberOfElements > 0)
	{
		std::vector<GeoLib::Point*>* points = new std::vector<GeoLib::Point*>();
		SHPObject* hSHPObject;

		for (int i = 0; i < numberOfElements; i++)
		{
			hSHPObject = SHPReadObject(hSHP,i);

			GeoLib::Point* pnt =
			        new GeoLib::Point( *(hSHPObject->padfX), *(hSHPObject->padfY),
			                           *(hSHPObject->padfZ) );
			points->push_back(pnt);
		}

		_geoObjects->addPointVec(points, listName);
		SHPDestroyObject(hSHPObject); // de-allocate SHPObject
	}
}

void SHPInterface::readStations(const SHPHandle &hSHP, int numberOfElements, std::string listName)
{
	if (numberOfElements > 0)
	{
		std::vector<GeoLib::Point*>* stations (new std::vector<GeoLib::Point*>);
		stations->reserve (numberOfElements);
		SHPObject* hSHPObject;

		for (int i = 0; i < numberOfElements; i++)
		{
			hSHPObject = SHPReadObject(hSHP,i);
			GeoLib::Station* stn =
			        GeoLib::Station::createStation( number2str(i), *(hSHPObject->padfX),
			                                        *(hSHPObject->padfY),
			                                        *(hSHPObject->padfZ) );
			stations->push_back(stn);
		}

		_geoObjects->addStationVec(stations, listName);
		SHPDestroyObject(hSHPObject); // de-allocate SHPObject
	}
}

void SHPInterface::readPolylines(const SHPHandle &hSHP, int numberOfElements, std::string listName)
{
	size_t noOfPoints = 0, noOfParts = 0;
	std::vector<GeoLib::Point*>* points = new std::vector<GeoLib::Point*>();
	std::vector<GeoLib::Polyline*>* lines = new std::vector<GeoLib::Polyline*>();
	SHPObject* hSHPObject;

	// for each polyline)
	for (int i = 0; i < numberOfElements; i++)
	{
		hSHPObject = SHPReadObject(hSHP,i);
		noOfPoints = hSHPObject->nVertices;
		noOfParts  = hSHPObject->nParts;

		for (size_t p = 0; p < noOfParts; p++)
		{
			int firstPnt = *(hSHPObject->panPartStart + p);
			int lastPnt  =
			        (p < (noOfParts - 1)) ? *(hSHPObject->panPartStart + p + 1) : noOfPoints;

			GeoLib::Polyline* line = new GeoLib::Polyline(*points);

			// for each point in that polyline
			for (int j = firstPnt; j < lastPnt; j++)
			{
				GeoLib::Point* pnt =
				        new GeoLib::Point( *(hSHPObject->padfX + j),
				                           *(hSHPObject->padfY + j),
				                           *(hSHPObject->padfZ + j) );
				points->push_back(pnt);
				line->addPoint(points->size() - 1);
			}

			// add polyline to polyline vector
			lines->push_back(line);
		}
	}

	if (numberOfElements > 0)
	{
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

	const std::vector<GeoLib::Polyline*>* polylines (_geoObjects->getPolylineVec(listName));
	std::vector<GeoLib::Surface*>* sfc_vec(new std::vector<GeoLib::Surface*>);

	for (std::vector<GeoLib::Polyline*>::const_iterator poly_it (polylines->begin());
	     poly_it != polylines->end(); poly_it++)
	{
		std::cout << "triangulation of Polygon with " << (*poly_it)->getNumberOfPoints() <<
		" points ... " << std::flush;
		sfc_vec->push_back(GeoLib::Surface::createSurface(*(*poly_it)));
	}

	if (!sfc_vec->empty())
		_geoObjects->addSurfaceVec(sfc_vec, listName);
}

void SHPInterface::adjustPolylines (std::vector<GeoLib::Polyline*>* lines,
                                    std::vector<size_t> id_map)

{
	for (size_t i = 0; i < lines->size(); i++)
	{
		GeoLib::Polyline* line( (*lines)[i] );
		size_t previous_pnt_id (std::numeric_limits<size_t>::max());

		for (size_t j = 0; j < line->getNumberOfPoints(); j++) {
			size_t jth_pnt_id(id_map[line->getPointID(j)]);
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
