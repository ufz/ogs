/**
 * \file
 * \author Thomas Fischer
 * \date   2010-12-01
 * \brief  Implementation of GEOTYPE enumeration helper functions.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GeoType.h"

#include <cstdlib>

namespace GeoLib
{

std::string convertGeoTypeToString (GEOTYPE geo_type)
{
    switch (geo_type)
    {
    case GEOTYPE::POINT:    return "POINT";
    case GEOTYPE::POLYLINE: return "POLYLINE";
    case GEOTYPE::SURFACE:  return "SURFACE";
    }

    std::abort(); // Cannot happen, because switch covers all cases.
                  // Used to silence compiler warning.
}

} // end namespace GeoLib
