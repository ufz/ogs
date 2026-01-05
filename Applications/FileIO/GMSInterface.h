// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

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
 * \brief Manages the import and export of Aquaveo GMS files into and out of
 * GeoLib.
 *
 * This class currently supports reading and writing ASCII borehole files as
 * well as (partially) reading mesh files. The 3dm-mesh-file-reader is based on
 * example meshes and does currently only support the following element types:
 * E4T (tetrahedra), E4P/E5P (pyramids) and E6W (wedges/prisms). Not supported
 * are E8H (Hex), E4Q (Quad), E3T (Tri) as well as higher order elements. Please
 * refer to the file format documentation of GMS for details.
 */
class GMSInterface final
{
public:
    /// Exports borehole data from all boreholes in a list to a file in
    /// GMS-format. (Note: there are some hardcoded tmp-files in the method that
    /// you might need to change!)
    static void writeBoreholesToGMS(const std::vector<GeoLib::Point*>* stations,
                                    const std::string& filename);

    /// Imports borehole data from a file in GMS-format.
    static int readBoreholesFromGMS(std::vector<GeoLib::Point*>& boreholes,
                                    const std::string& filename);

    /// Reads a GMS *.3dm file and converts it to an CFEMesh.
    static MeshLib::Mesh* readMesh(const std::string& filename);
};

}  // namespace FileIO
