// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "GeoType.h"

#include <cstdlib>

#include "BaseLib/Error.h"

namespace GeoLib
{
std::string convertGeoTypeToString(GEOTYPE geo_type)
{
    switch (geo_type)
    {
        case GEOTYPE::POINT:
            return "POINT";
        case GEOTYPE::POLYLINE:
            return "POLYLINE";
        case GEOTYPE::SURFACE:
            return "SURFACE";
    }

    // Cannot happen, because switch covers all cases.
    // Used to silence compiler warning.
    OGS_FATAL("convertGeoTypeToString(): Given geo type is not supported");
}

}  // end namespace GeoLib
