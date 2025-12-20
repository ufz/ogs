// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string>

namespace GeoLib
{

enum class GEOTYPE
{
    POINT,
    POLYLINE,
    SURFACE
};

std::string convertGeoTypeToString(GEOTYPE geo_type);

}  // end namespace GeoLib
