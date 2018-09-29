/**
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>

namespace GeoLib
{
    class GEOObjects;
}

namespace FileIO
{
    void
    readGeometryFromFile(std::string const& fname, GeoLib::GEOObjects & geo_objs);
}
