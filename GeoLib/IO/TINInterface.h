/**
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <string>
#include <vector>

#include "GeoLib/Point.h"
#include "GeoLib/PointVec.h"

namespace GeoLib
{
class Surface;

namespace IO
{

/**
 * Interface for reading and writing Triangulated Irregular Network (TIN) file formats.
 */
class TINInterface
{
public:
    /**
     * Reads TIN file
     * @param fname    TIN file name
     * @param pnt_vec  a point vector to which triangle vertexes are added
     * @param errors   a vector of error messages
     * @return a pointer to a GeoLib::Surface object created from TIN data. nullptr is returned if it fails to read the file.
     */
    static GeoLib::Surface* readTIN(std::string const& fname,
                                    GeoLib::PointVec &pnt_vec,
                                    std::vector<std::string>* errors = nullptr);

    /**
     * Writes surface data into TIN file
     * @param surface    surface object
     * @param file_name  TIN file name
     */
    static void writeSurfaceAsTIN(GeoLib::Surface const& surface, std::string const& file_name);
};

} // end namespace IO
} // end namespace GeoLib
