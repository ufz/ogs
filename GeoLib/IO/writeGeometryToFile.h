/**
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef WRITEGEOMETRYTOFILE_H
#define WRITEGEOMETRYTOFILE_H

#include <string>

namespace GeoLib
{
class GEOObjects;

namespace IO
{

/// Write geometry given by the \c geo_objs object and specified by the name
/// stored in param \c geo_name either to a gml or a gli file. If the extension
/// given in the \c fname parameter is "gml" or "GML" a gml file is written. In
/// case the extension is "gli" or "GLI" a gli file is written.
void writeGeometryToFile(std::string const& geo_name,
    GeoLib::GEOObjects& geo_objs, std::string const& fname);
}
}

#endif // WRITEGEOMETRYTOFILE_H
