/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-17
 * \brief  Definition of the GEOTYPE enumeration.
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

namespace GeoLib {

/**
 * \ingroup GeoLib
 */

enum class GEOTYPE {
    POINT,
    POLYLINE,
    SURFACE
};

std::string convertGeoTypeToString (GEOTYPE geo_type);

} // end namespace GeoLib
