/**
 * \file
 * \author Karsten Rink
 * \date   2010-06-08
 * \brief  Definition of the GMSInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * @file GMSInterface.h
 * @date 2010-06-08
 * @author Lars Bilke
 */

#ifndef GMSINTERFACE_H_
#define GMSINTERFACE_H_

#include <list>
#include <string>
#include <vector>

namespace GeoLib
{
class Station;
class StationBorehole;
class Point;
}

namespace MeshLib
{
class Mesh;
}

namespace FileIO
{

/**
 * \brief Manages the import and export of Aquaveo GMS files into and out of GeoLib.
 *
 * This class currently supports reading and writing ASCII borehole files as well as
 * (partially) reading mesh files.
 * The 3dm-mesh-file-reader is based on example meshes and does currently only support
 * the following element types: E4T (tetrahedra), E4P/E5P (pyramids) and E6W (wedges/prisms).
 * Not supported are E8H (Hex), E4Q (Quad), E3T (Tri) as well as higher order elements.
 * Please refer to the file format documentation of GMS for details.
 */
class GMSInterface
{
public:
    /// Exports borehole data from all boreholes in a list to a file in GMS-format. (Note: there are some hardcoded tmp-files in the method that you might need to change!)
    static void writeBoreholesToGMS(const std::vector<GeoLib::Point*>* stations,
                                    const std::string &filename);

    /// Imports borehole data from a file in GMS-format.
    static int readBoreholesFromGMS(std::vector<GeoLib::Point*>* boreholes,
                                    const std::string &filename);

    /// Writes a file that assigns each soilID-index in the GMS export file a name.
    static int writeSoilIDTable(const std::vector<std::string> &soilID,
                                const std::string &filename);

    /// Reads a GMS *.3dm file and converts it to an CFEMesh.
    static MeshLib::Mesh* readGMS3DMMesh(const std::string &file_name);

private:
    /**
     * \brief Reads SoilIDs for Borehole export from an external file
     *
     * The method expects a file with the name of one stratigraphic layer at each line. These layers are assigned
     * ascending IDs, i.e. the first name gets index 0, the second line gets index 1, etc.
     * \return An array with the names of the stratigraphic layers in which the index for each string equals its ID.
     */
    static std::vector<std::string> readSoilIDfromFile(const std::string &filename);

    /// Finds the ID assigned to soilName or creates a new one ( this method is called from writeBoreholeToGMS() )
    static std::size_t getSoilID(std::vector<std::string> &soilID, std::string &soilName);
};

}

#endif /* GMSINTERFACE_H_ */
