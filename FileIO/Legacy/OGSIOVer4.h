/**
 * \file
 * \author Thomas Fischer
 * \date   2010-01-14
 * \brief  Definition of the OGSIOVer4 class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OGSIOVER4_H_
#define OGSIOVER4_H_

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "Point.h"
#include "Polyline.h"
#include "StringTools.h"
#include "Surface.h"

// forward declaration
namespace GeoLib
{
class GEOObjects;
}

namespace FileIO
{
namespace Legacy {

/** Interface for handling geometry from OGS-5 and below (*.gli files) */

/** Reads geometric objects from file in gli format */
bool readGLIFileV4 (const std::string& fname, GeoLib::GEOObjects* geo, std::string& unique_name, std::vector<std::string>& errors);

/** Writes geometric objects from a specific geometry to a gli-file */
void writeGLIFileV4 (const std::string& fname,
                     const std::string& proj_name,
                     const GeoLib::GEOObjects& geo);

/** Writes all geometric information to a gli-file */
void writeAllDataToGLIFileV4 (const std::string& fname, const GeoLib::GEOObjects& geo);

}
} // end namespace FileIO

#endif /* OGSIOVER4_H_ */
