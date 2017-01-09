/**
 * \file
 * \author Karsten Rink
 * \date   2010-01-25
 * \brief  Implementation of the SHPInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * @file SHPInterface.h
 * @date 25/01/2010
 * @author Karsten Rink
 */

#ifndef SHPINTERFACE_H
#define SHPINTERFACE_H

#include <string>
#include <vector>

//ShapeLib includes
#include <shapefil.h>


namespace GeoLib {
    class GEOObjects;
    class Polyline;
}

namespace MeshLib {
    class Mesh;
}


namespace FileIO
{

/**
 * \brief Manages the import of ESRI shape files into GeoLib.
 */
class SHPInterface
{
public:
    /// Connection between ESRI type system for shape files and OGS GeoLib.
    enum class OGSType
    {
        UNDEFINED   = 0,
        POINT       = 1,
        STATION     = 2,
        POLYLINE    = 3,
        POLYGON     = 4
    };

    /// Constructor
    SHPInterface(GeoLib::GEOObjects& geoObjects) : _geoObjects(geoObjects) {}

    /// Reads the header of the shape file.
    bool readSHPInfo(const std::string &filename, int &shapeType, int &numberOfEntities);

    /// Reads data from the shape file.
    void readSHPFile(const std::string &filename, OGSType choice, const std::string &listName);

    /// Writes a 2D mesh into a shapefile using one polygon for every element
    /// (based on request by AS, open for discussion)
    static bool write2dMeshToSHP(const std::string &file_name, const MeshLib::Mesh &mesh);


private:
    /// Reads points into a vector of Point objects.
    void readPoints    (const SHPHandle &hSHP, int numberOfElements, std::string listName);

    /// Reads points into a vector of Point objects and marks them as Station.
    void readStations  (const SHPHandle &hSHP, int numberOfElements, std::string listName);

    /// Reads lines into a vector of Polyline objects.
    void readPolylines (const SHPHandle &hSHP, int numberOfElements, std::string listName);

    /// Reads lines into a vector of Polyline and Surface objects.
    void readPolygons  (const SHPHandle &hSHP, int numberOfElements, const std::string& listName);

    GeoLib::GEOObjects& _geoObjects;
};

}

#endif //SHPINTERFACE_H
